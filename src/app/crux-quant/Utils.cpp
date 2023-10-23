#pragma once

#include "Utils.h"

#include <algorithm>
#include <cfloat>  // For DBL_MAX
#include <cmath>
#include <iostream>
#include <map>
#include <mutex>  // Include the mutex header
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <tuple>



#include "csv.h"
#include "correlation.h"


using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

std::mutex mtx;  // Declare a mutex

const int MaxThreads = 4;

typedef pwiz::msdata::SpectrumListPtr SpectrumListPtr;

namespace CruxQuant {

SpectrumListPtr loadSpectra(const string& file, int msLevel) {
    try {
        MSDataFile msd(file);
        SpectrumListPtr originalSpectrumList = msd.run.spectrumListPtr;
        if (!originalSpectrumList) {
            carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
        }
        SpectrumListSimplePtr filteredSpectrumList(new SpectrumListSimple);

        for (size_t i = 0; i < originalSpectrumList->size(); ++i) {
            SpectrumPtr spectrum = originalSpectrumList->spectrum(i);

            int spectrumMSLevel = spectrum->cvParam(MS_ms_level).valueAs<int>();
            if (spectrumMSLevel == msLevel) {
                // Add the spectrum to the filtered list
                filteredSpectrumList->spectra.push_back(spectrum);
            }
        }

        return filteredSpectrumList;
    } catch (exception& e) {
        carp(CARP_INFO, "Error:  %s", e.what());
        return nullptr;
    }
}

IndexedSpectralResults indexedMassSpectralPeaks(SpectrumListPtr spectrum_collection, const string& spectra_file) {
    string _spectra_file(spectra_file);

    map<int, map<int, IndexedMassSpectralPeak>> _indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;

    _ms1Scans[_spectra_file] = vector<Ms1ScanInfo>();

    IndexedSpectralResults index_results{_indexedPeaks, _ms1Scans};

    if (!spectrum_collection) {
        return index_results;
    }

    int _scanIndex = 0;
    int _oneBasedScanNumber = 1;

    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        SpectrumPtr spectrum = spectrum_collection->spectrum(i);

        if (spectrum) {
            BinaryDataArrayPtr mzs = spectrum->getMZArray();
            BinaryDataArrayPtr intensities = spectrum->getIntensityArray();

            int scanIndex;
            int oneBasedScanNumber;
            std::string scanId = getScanID(spectrum->id);
            if (scanId.empty()) {
                scanIndex = _scanIndex;
                oneBasedScanNumber = _oneBasedScanNumber;
            } else {
                scanIndex = std::stoi(scanId);
                oneBasedScanNumber = scanIndex + 1;
            }

            double retentionTime = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

            if (mzs && intensities) {
                const std::vector<double>& mzArray = mzs->data;
                const std::vector<double>& intensityArray = intensities->data;
                for (size_t j = 0; j < mzArray.size(); ++j) {
                    FLOAT_T mz = mzArray[j];
                    int roundedMz = static_cast<int>(std::round(mz * BINS_PER_DALTON));
                    IndexedMassSpectralPeak spec_data(
                        mz,                 // mz value
                        intensityArray[j],  // intensity
                        scanIndex,          // zeroBasedMs1ScanIndex
                        retentionTime);

                    auto& indexedPeaks = index_results._indexedPeaks;
                    auto it = indexedPeaks.find(roundedMz);
                    if (it == indexedPeaks.end()) {
                        map<int, IndexedMassSpectralPeak> tmp;
                        tmp.insert({scanIndex, spec_data});
                        indexedPeaks[roundedMz] = tmp;
                    } else {
                        it->second.insert({scanIndex, spec_data});
                    }

                    Ms1ScanInfo scan = {oneBasedScanNumber, scanIndex, retentionTime};
                    index_results._ms1Scans[spectra_file].push_back(scan);
                }
            }

            _scanIndex++;
            _oneBasedScanNumber++;
        }
    }
    return index_results;
}

map<int, PSM> create_psm_map(const string& psm_file) {
    carp(CARP_INFO, "loading psm data ...");
    map<int, PSM> psm_datum;

    io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>> matchFileReader(psm_file);

    string sequence_col;
    int scan_col, charge_col;
    double peptide_mass_col, spectrum_precursor_mz_col;

    matchFileReader.read_header(io::ignore_extra_column, "scan", "charge", "spectrum precursor m/z", "peptide mass", "sequence");

    while (matchFileReader.read_row(scan_col, charge_col, spectrum_precursor_mz_col, peptide_mass_col, sequence_col)) {
        PSM psm = {
            sequence_col,
            scan_col,
            charge_col,
            peptide_mass_col,
            spectrum_precursor_mz_col};
        psm_datum[scan_col] = psm;
    }

    return psm_datum;
}
// Make this a multithreaded process
vector<Identification> createIdentifications(const map<int, PSM>& psm_datum, const string& spectra_file, SpectrumListPtr spectrum_collection) {
    carp(CARP_INFO, "Creating indentifications, this may take a bit of time, do not terminate the process...");

    vector<Identification> allIdentifications;
    string _spectra_file(spectra_file);

    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        SpectrumPtr spectrum = spectrum_collection->spectrum(i);
        if (spectrum) {
            int scanIndex;
            std::string scanId = getScanID(spectrum->id);
            if (scanId.empty()) {
                continue;
            } else {
                scanIndex = std::stoi(scanId);
            }

            auto it = psm_datum.find(scanIndex);

            if (it != psm_datum.end()) {
                double retentionTimeInSeconds = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

                FLOAT_T retentionTimeInMinutes = retentionTimeInSeconds / 60.0;

                Identification identification;

                identification.sequence = it->second.sequence_col;
                identification.monoIsotopicMass = it->second.peptide_mass_col;
                identification.charge = it->second.charge_col;
                identification.peptideMass = it->second.peptide_mass_col;
                identification.precursorCharge = it->second.spectrum_precursor_mz_col;
                identification.spectralFile = _spectra_file;
                identification.ms2RetentionTimeInMinutes = retentionTimeInMinutes;
                identification.scanId = it->second.scan_col;
                allIdentifications.push_back(identification);
            }
        }
    }

    return allIdentifications;
}

string calcFormula(string seq) {
    int H = 2;
    int C = 0;
    int N = 0;
    int O = 1;
    int S = 0;
    char str[128];
    string s;
    string mod;
    int x;

    size_t i;
    for (i = 0; i < seq.size(); i++) {
        switch (seq[i]) {
            case 'A':
                C += 3;
                H += 5;
                N += 1;
                O += 1;
                break;
            case 'R':
                C += 6;
                H += 12;
                N += 4;
                O += 1;
                break;
            case 'N':
                C += 4;
                H += 6;
                N += 2;
                O += 2;
                break;
            case 'D':
                C += 4;
                H += 5;
                N += 1;
                O += 3;
                break;
            case 'C':
                C += 5;
                H += 8;
                N += 2;
                O += 2;
                S += 1;
                break;  // carbamidomethylation of C
            case 'Q':
                C += 5;
                H += 8;
                N += 2;
                O += 2;
                break;
            case 'E':
                C += 5;
                H += 7;
                N += 1;
                O += 3;
                break;
            case 'G':
                C += 2;
                H += 3;
                N += 1;
                O += 1;
                break;
            case 'H':
                C += 6;
                H += 7;
                N += 3;
                O += 1;
                break;
            case 'I':
            case 'L':
                C += 6;
                H += 11;
                N += 1;
                O += 1;
                break;
            case 'K':
                C += 6;
                H += 12;
                N += 2;
                O += 1;
                break;
            case 'M':
                C += 5;
                H += 9;
                N += 1;
                O += 1;
                S += 1;
                break;
            case 'F':
                C += 9;
                H += 9;
                N += 1;
                O += 1;
                break;
            case 'P':
                C += 5;
                H += 7;
                N += 1;
                O += 1;
                break;
            case 'S':
                C += 3;
                H += 5;
                N += 1;
                O += 2;
                break;
            case 'T':
                C += 4;
                H += 7;
                N += 1;
                O += 2;
                break;
            case 'W':
                C += 11;
                H += 10;
                N += 2;
                O += 1;
                break;
            case 'Y':
                C += 9;
                H += 9;
                N += 1;
                O += 2;
                break;
            case 'V':
                C += 5;
                H += 9;
                N += 1;
                O += 1;
                break;
            case '[':
                mod.clear();
                break;
            case '1':
                mod += '1';
                break;
            case '2':
                mod += '2';
                break;
            case '3':
                mod += '3';
                break;
            case '4':
                mod += '4';
                break;
            case '5':
                mod += '5';
                break;
            case '6':
                mod += '6';
                break;
            case '7':
                mod += '7';
                break;
            case '8':
                mod += '8';
                break;
            case '9':
                mod += '9';
                break;
            case ']':
                x = atoi(mod.c_str());
                if (x == 147)
                    O += 1;
                break;
            default:
                break;
        }
    }

    if (S > 0) {
        sprintf(str, "C%dH%dN%dO%dS%d", C, H, N, O, S);
    } else {
        sprintf(str, "C%dH%dN%dO%d", C, H, N, O);
    }

    s = str;
    return s;
}

unordered_map<string, vector<pair<double, double>>> calculateTheoreticalIsotopeDistributions(const vector<Identification>& allIdentifications) {
    unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution;

    for (const auto& identification : allIdentifications) {
        string peptide_sequence = identification.sequence;
        int charge = identification.charge;
        double peptide_mass = identification.peptideMass;

        if (modifiedSequenceToIsotopicDistribution.find(peptide_sequence) != modifiedSequenceToIsotopicDistribution.end()) {
            continue;
        }

        vector<pair<double, double>> isotopicMassesAndNormalizedAbundances;

        string formula = calcFormula(peptide_sequence);
        char* char_array = new char[formula.length() + 1];
        char* fn = nullptr;
        strcpy(char_array, formula.c_str());
        CMercury8 dist(fn);
        // dist.Echo(true);
        dist.GoMercury(char_array, charge);
        vector<double> masses;
        vector<double> abundances;
        for (auto i : dist.FixedData) {
            masses.push_back(i.mass);
            abundances.push_back(i.data);
        }

        double highestAbundance = *std::max_element(abundances.begin(), abundances.end());

        for (int i = 0; i < masses.size(); i++) {
            // Calculate the expected isotopic mass shift for this peptide
            masses[i] -= peptide_mass;

            // Normalize the abundance of each isotope
            abundances[i] /= highestAbundance;

            // Look for these isotopes
            if (isotopicMassesAndNormalizedAbundances.size() < NUM_ISOTOPES_REQUIRED || abundances[i] > 0.1) {
                isotopicMassesAndNormalizedAbundances.push_back(std::make_pair(masses[i], abundances[i]));
            }
        }
        modifiedSequenceToIsotopicDistribution[peptide_sequence] = isotopicMassesAndNormalizedAbundances;
        delete[] char_array;
    }

    return modifiedSequenceToIsotopicDistribution;
}

void setPeakFindingMass(vector<Identification>& allIdentifications, unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution) {
    for (auto& identification : allIdentifications) {
        const string& sequence = identification.sequence;

        // Find the isotope where normalized abundance is 1
        double mostAbundantIsotopeShift = 0.0;

        const auto& isotopicDistribution = modifiedSequenceToIsotopicDistribution.find(sequence);
        if (isotopicDistribution != modifiedSequenceToIsotopicDistribution.end()) {
            for (const auto& item : isotopicDistribution->second) {
                if (item.second == 1.0) {
                    mostAbundantIsotopeShift = item.first;
                    break;
                }
            }
        }

        identification.peakFindingMass = identification.monoIsotopicMass + mostAbundantIsotopeShift;
    }
}

vector<double> createChargeStates(const vector<Identification>& allIdentifications) {
    auto minChargeState = std::numeric_limits<double>::max();
    auto maxChargeState = std::numeric_limits<double>::lowest();
    for (const auto& identification : allIdentifications) {
        minChargeState = std::min(minChargeState, identification.precursorCharge);
        maxChargeState = std::max(maxChargeState, identification.precursorCharge);
    }
    vector<double> chargeStates;

    for (double charge = minChargeState; charge <= maxChargeState; charge++) {
        chargeStates.push_back(charge);
    }

    return chargeStates;
}

void processRange(int start, int end,
                  const vector<Identification>& ms2IdsForThisFile,
                  const string& spectralFile,
                  const vector<double>& chargeStates,
                  vector<ChromatographicPeak>& chromatographicPeaks,
                  PpmTolerance& peakfindingTol,
                  unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
                  map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
                  PpmTolerance& ppmTolerance,
                  unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution) {
    for (int i = start; i < end; ++i) {
        const Identification& identification = ms2IdsForThisFile[i];
        {
            std::lock_guard<std::mutex> lock(mtx);
            ChromatographicPeak msmsFeature(identification, false, spectralFile);

            msmsFeature._index = i;
            chromatographicPeaks.push_back(msmsFeature);
            for (const auto& chargeState : chargeStates) {
                if (ID_SPECIFIC_CHARGE_STATE && chargeState != identification.precursorCharge) {
                    continue;
                }

                vector<IndexedMassSpectralPeak*> xic = peakFind(
                    identification.ms2RetentionTimeInMinutes,  // This value may not be in minutes, check and convert
                    identification.peakFindingMass,
                    chargeState,
                    identification.spectralFile,
                    peakfindingTol,
                    _ms1Scans,
                    indexedPeaks);

                // Sort the xic vector based on the RetentionTime
                std::sort(xic.begin(), xic.end(), [](const IndexedMassSpectralPeak* a, const IndexedMassSpectralPeak* b) {
                    return a->retentionTime < b->retentionTime;
                });

                xic.erase(std::remove_if(xic.begin(), xic.end(), [&](const IndexedMassSpectralPeak* p) {
                              return !ppmTolerance.Within(toMass(p->mz, chargeState), identification.peakFindingMass);
                          }),
                          xic.end());

                vector<IsotopicEnvelope> isotopicEnvelopes = getIsotopicEnvelopes(xic, identification, chargeState, modifiedSequenceToIsotopicDistribution, indexedPeaks);

                msmsFeature.isotopicEnvelopes.insert(msmsFeature.isotopicEnvelopes.end(), isotopicEnvelopes.begin(), isotopicEnvelopes.end());
            }

            msmsFeature.calculateIntensityForThisFeature(INTEGRATE);

            // TO DO - Continue from here

        }
    }
}

void quantifyMs2IdentifiedPeptides(
    string spectraFile,
    const vector<Identification>& allIdentifications,
    const vector<double>& chargeStates,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution) {
    carp(CARP_INFO, "Quantifying MS2, this may take some time...");
    vector<Identification> ms2IdsForThisFile;

    // Use std::copy_if to filter the identifications
    std::copy_if(
        allIdentifications.begin(),
        allIdentifications.end(),
        std::back_inserter(ms2IdsForThisFile),
        [&spectraFile](const Identification& id) {
            return id.spectralFile == spectraFile;
        });

    // Check if ms2IdsForThisFile is empty
    if (ms2IdsForThisFile.empty()) {
        return;  // Early return
    }

    PpmTolerance peakfindingTol(PEAK_FINDING_PPM_TOLERANCE);  // Peak finding tolerance is generally higher than ppmTolerance
    PpmTolerance ppmTolerance(PPM_TOLERANCE);
    int totalCount = ms2IdsForThisFile.size();

    vector<ChromatographicPeak> chromatographicPeaks;

    vector<std::thread> threads;

    int chunkSize = totalCount / MaxThreads;
    int remainder = totalCount % MaxThreads;

    // Create and start the threads
    for (int i = 0; i < MaxThreads; ++i) {
        int start = i * chunkSize;
        int end = (i == MaxThreads - 1) ? (start + chunkSize + remainder) : (start + chunkSize);

        threads.push_back(std::thread([start, end,
                                       &ms2IdsForThisFile,
                                       &spectraFile,
                                       &chargeStates,
                                       &chromatographicPeaks,
                                       &peakfindingTol,
                                       &_ms1Scans,
                                       &indexedPeaks,
                                       &ppmTolerance,
                                       &modifiedSequenceToIsotopicDistribution]() {
            processRange(start, end,
                         ms2IdsForThisFile,
                         spectraFile,
                         chargeStates,
                         chromatographicPeaks,
                         peakfindingTol,
                         _ms1Scans,
                         indexedPeaks,
                         ppmTolerance,
                         modifiedSequenceToIsotopicDistribution);
        }));
    }

    // Join the threads
    for (auto& thread : threads) {
        thread.join();
    }
}

double toMz(double mass, int charge) {
    return mass / std::abs(charge) + std::copysign(PROTONMASS, charge);
}

double toMass(double massToChargeRatio, int charge) {
    return std::abs(charge) * massToChargeRatio - charge * PROTONMASS;
}

IndexedMassSpectralPeak* getIndexedPeak(double theorMass, int zeroBasedScanIndex, PpmTolerance tolerance, int chargeState, map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
    IndexedMassSpectralPeak* bestPeak = nullptr;

    // Calculate the maximum value using tolerance
    double maxValue = tolerance.GetMaximumValue(theorMass);
    double maxMzValue = toMz(maxValue, chargeState);
    // Calculate ceilingMz by rounding maxMzValue up to the nearest integer
    int ceilingMz = static_cast<int>(std::ceil(maxMzValue * BINS_PER_DALTON));

    double minValue = tolerance.GetMinimumValue(theorMass);
    double minMzValue = toMz(minValue, chargeState);
    int floorMz = static_cast<int>(std::floor(minMzValue * BINS_PER_DALTON));

    for (int j = floorMz; j <= ceilingMz; j++) {
        map<int, IndexedMassSpectralPeak> bin = indexedPeaks[j];

        auto startIterator = bin.find(zeroBasedScanIndex);

        for (auto it = startIterator; it != bin.end(); ++it) {
            auto _peak = it->second;
            if (_peak.zeroBasedMs1ScanIndex > zeroBasedScanIndex) {
                break;
            }

            double expMass = toMass(_peak.mz, chargeState);

            if (tolerance.Within(expMass, theorMass) && _peak.zeroBasedMs1ScanIndex == zeroBasedScanIndex && (bestPeak == nullptr || std::abs(expMass - theorMass) < std::abs(toMass(bestPeak->mz, chargeState) - theorMass))) {
                bestPeak = &_peak;
            }
        }
    }
    return bestPeak;
}
void ChromatographicPeak::calculateIntensityForThisFeature(bool integrate) {
    if (!isotopicEnvelopes.empty()) {
        // Find the IsotopicEnvelope with the maximum intensity (apex)
        auto maxIntensityEnvelopes = std::max_element(
            isotopicEnvelopes.begin(), isotopicEnvelopes.end(),
            [](const IsotopicEnvelope& a, const IsotopicEnvelope& b) {
                return a.intensity < b.intensity;
            });

        if (maxIntensityEnvelopes != isotopicEnvelopes.end()) {
            apex = *maxIntensityEnvelopes;

            if (integrate) {
                // Calculate intensity by summing up all IsotopicEnvelope intensities
                intensity = std::accumulate(
                    isotopicEnvelopes.begin(), isotopicEnvelopes.end(), 0.0,
                    [](double sum, const IsotopicEnvelope& envelope) {
                        return sum + envelope.intensity;
                    });
            } else {
                // Set intensity to the apex intensity
                intensity = apex.intensity;
            }

            // Calculate massError for each Identification
            massError = std::numeric_limits<double>::quiet_NaN();  // Set massError to NaN initially

            for (const Identification& id : identifications) {
                double massErrorForId = ((toMass(apex.indexedPeak.mz, apex.chargeState) - id.peakFindingMass) / id.peakFindingMass) * 1e6;

                if (std::isnan(massError) || std::abs(massErrorForId) < std::abs(massError)) {
                    massError = massErrorForId;
                }
            }
            std::unordered_set<int> uniqueChargeStates;
            for (const IsotopicEnvelope& envelope : isotopicEnvelopes) {
                uniqueChargeStates.insert(envelope.chargeState);
            }
            numChargeStatesObserved = static_cast<int>(uniqueChargeStates.size());
        }
    } else {
        // No isotopicEnvelopes, so set intensity to 0 and massError to NaN
        intensity = 0.0;
        massError = std::numeric_limits<double>::quiet_NaN();
        numChargeStatesObserved = 0;
        apex = {};  // Reset apex to default-constructed IsotopicEnvelope
    }
}

vector<IndexedMassSpectralPeak*> peakFind(double idRetentionTime,
                                          double mass,
                                          int charge,
                                          const string& spectra_file,
                                          PpmTolerance tolerance,
                                          unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
                                          map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
    vector<IndexedMassSpectralPeak*> xic;
    vector<Ms1ScanInfo> ms1Scans = _ms1Scans[spectra_file];
    int precursorScanIndex = -1;
    for (const Ms1ScanInfo& ms1Scan : ms1Scans) {
        if (ms1Scan.retentionTime < idRetentionTime) {
            precursorScanIndex = ms1Scan.zeroBasedMs1ScanIndex;
        } else {
            break;
        }
    }

    // go right
    int missedScans = 0;
    for (int t = precursorScanIndex; t < ms1Scans.size(); t++) {
        auto peak = getIndexedPeak(mass, t, tolerance, charge, indexedPeaks);
        if (peak == nullptr && t != precursorScanIndex) {
            missedScans++;
        } else if (peak != nullptr) {
            missedScans = 0;
            xic.push_back(peak);
        }

        if (missedScans > MISSED_SCANS_ALLOWED) {
            break;
        }
    }

    // go left
    missedScans = 0;
    for (int t = precursorScanIndex - 1; t >= 0; t--) {
        auto peak = getIndexedPeak(mass, t, tolerance, charge, indexedPeaks);

        if (peak == nullptr && t != precursorScanIndex) {
            missedScans++;
        } else if (peak != nullptr) {
            missedScans = 0;
            xic.push_back(peak);
        }

        if (missedScans > MISSED_SCANS_ALLOWED) {
            break;
        }
    }

    std::sort(xic.begin(), xic.end(), [](const IndexedMassSpectralPeak* x, const IndexedMassSpectralPeak* y) {
        return x->retentionTime < y->retentionTime;
    });

    return xic;
}

string getScanID(string spectrum_id) {
    size_t scanPos = spectrum_id.find("scan=");
    if (scanPos != std::string::npos) {
        return spectrum_id.substr(scanPos + 5);
    }
    return "";
}

vector<IsotopicEnvelope> getIsotopicEnvelopes(
    const vector<IndexedMassSpectralPeak*>& xic,
    const Identification& identification,
    const int chargeState,
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
    vector<IsotopicEnvelope> isotopicEnvelopes;
    vector<pair<double, double>> isotopeMassShifts = modifiedSequenceToIsotopicDistribution[identification.sequence];

    if (isotopeMassShifts.size() < NUM_ISOTOPES_REQUIRED) {
        return isotopicEnvelopes;
    }

    PpmTolerance isotopeTolerance(ISOTOPE_TOLERANCE_PPM);

    vector<double> experimentalIsotopeIntensities(isotopeMassShifts.size());

    vector<double> theoreticalIsotopeMassShifts;
    vector<double> theoreticalIsotopeAbundances;
    for (const auto& p : isotopeMassShifts) {
        theoreticalIsotopeMassShifts.push_back(p.first);
        theoreticalIsotopeAbundances.push_back(p.second);
    }

    int peakfindingMassIndex = static_cast<int>(round(identification.peakFindingMass - identification.monoIsotopicMass));

    // For each peak in the XIC, we consider the possibility that there was an off-by-one or missed monoisotopic mass
    // error in peak assignment / deconvolution. The -1 key in this dictionary corresponds to a negative off-by-one error, the
    // +1 key corresponds to a positive off-by-one error, and the 0 key corresponds to accurate assignment/deconvolution.
    map<int, vector<IsotopePeak>> massShiftToIsotopePeaks = {
        {-1, {}},
        {0, {}},
        {1, {}}};

    vector<int> directions = {-1, 1};

    for (const IndexedMassSpectralPeak* peak : xic) {
        std::fill_n(experimentalIsotopeIntensities.begin(), experimentalIsotopeIntensities.size(), 0);

        for (auto& kvp : massShiftToIsotopePeaks) {
            kvp.second.clear();
        }

        // isotope masses are calculated relative to the observed peak
        double observedMass = toMass(peak->mz, chargeState);
        double observedMassError = observedMass - identification.peakFindingMass;

        for (auto& shift : massShiftToIsotopePeaks) {
            // look for each isotope peak in the data
            // This is done by starting with the first isotope with mass less than the
            // peak finding (most abundant) mass, then working backwards to find every isotope
            // with mass < most abundant mass. Once an expected isotopic peak can not be found,
            // the loop breaks and we begin working our way forward, starting with the peak finding
            // mass and locating every peak with mass > peak finding mass. Once an expected isotopic
            // peak can not be found, it is assumed that we have located every isotope present and the loop breaks.

            for (int direction : directions) {
                int start = (direction == -1) ? peakfindingMassIndex - 1 : peakfindingMassIndex;

                for (int i = start; i < theoreticalIsotopeAbundances.size() && i >= 0; i += direction) {
                    double isotopeMass = identification.monoIsotopicMass + observedMassError +
                                         theoreticalIsotopeMassShifts[i] + shift.first * C13MinusC12;
                    double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] * peak->intensity;

                    IndexedMassSpectralPeak* isotopePeak = getIndexedPeak(isotopeMass,
                                                                          peak->zeroBasedMs1ScanIndex, isotopeTolerance, chargeState, indexedPeaks);

                    if (isotopePeak == nullptr || isotopePeak->intensity < theoreticalIsotopeIntensity / 4.0 || isotopePeak->intensity > theoreticalIsotopeIntensity * 4.0) {
                        break;
                    }

                    shift.second.push_back(std::make_tuple(isotopePeak->intensity, theoreticalIsotopeIntensity, isotopeMass));
                    if (shift.first == 0) {
                        experimentalIsotopeIntensities[i] = isotopePeak->intensity;
                    }
                }
            }

            // check number of isotope peaks observed
            if (massShiftToIsotopePeaks[0].size() < NUM_ISOTOPES_REQUIRED) {
                continue;
            }

            // Check that the experimental envelope matches the theoretical
            if (checkIsotopicEnvelopeCorrelation(massShiftToIsotopePeaks, peak, chargeState, isotopeTolerance, indexedPeaks)) {
                for (size_t i = 0; i < experimentalIsotopeIntensities.size(); ++i) {
                    if (experimentalIsotopeIntensities[i] == 0) {
                        experimentalIsotopeIntensities[i] = theoreticalIsotopeAbundances[i] * experimentalIsotopeIntensities[peakfindingMassIndex];
                    }
                }
            }

            double sumExperimentalIntensities = std::accumulate(experimentalIsotopeIntensities.begin(), experimentalIsotopeIntensities.end(), 0.0);

            IsotopicEnvelope isotopicEnvelope(*peak, chargeState, sumExperimentalIntensities);

            isotopicEnvelopes.push_back(isotopicEnvelope);
        }
    }

    return isotopicEnvelopes;
}

filterResults filterMassShiftToIsotopePeaks(vector<IsotopePeak>& isotopePeaks) {
    std::vector<double> expIntensity;
    std::vector<double> theorIntensity;

    for (const IsotopePeak& peak : isotopePeaks) {
        double _expIntensity = std::get<0>(peak);
        expIntensity.push_back(_expIntensity);
        double _theorIntensity = std::get<1>(peak);
        theorIntensity.push_back(_theorIntensity);
    }

    filterResults results {
        expIntensity,
        theorIntensity
    };

    return results;
}

void setToNegativeOneIfNaN(double& value) {
    if (std::isnan(value)) {
        value = -1;
    }
}

bool checkIsotopicEnvelopeCorrelation(
    map<int, vector<IsotopePeak>>& massShiftToIsotopePeaks,
    const IndexedMassSpectralPeak* peak,
    int chargeState,
    PpmTolerance& isotopeTolerance,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks
    ) {
    filterResults results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[0]);
    double corr = Pearson(results.expIntensity, results.theorIntensity);

    for (auto& shift : massShiftToIsotopePeaks) {
        if (shift.second.empty()) {
            continue;
        }

        double unexpectedMass = DBL_MAX;  // Initialize with a large value
        for (const IsotopePeak& _peak : shift.second) {
            double theorMass = std::get<2>(_peak);
            unexpectedMass = std::min(unexpectedMass, theorMass);
        }
        unexpectedMass -= C13MinusC12;
        IndexedMassSpectralPeak* unexpectedPeak = getIndexedPeak(unexpectedMass, peak->zeroBasedMs1ScanIndex, isotopeTolerance, chargeState, indexedPeaks);

        if (unexpectedPeak == nullptr) {
            shift.second.push_back(std::make_tuple(0.0, 0.0, unexpectedMass));
        } else {
            shift.second.push_back(std::make_tuple(unexpectedPeak->intensity, 0.0, unexpectedMass));
        }
    }

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[0]);
    double corrWithPadding = Pearson(results.expIntensity, results.theorIntensity);

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[-1]);
    double corrShiftedLeft = Pearson(results.expIntensity, results.theorIntensity);

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[1]);
    double corrShiftedRight = Pearson(results.expIntensity, results.theorIntensity);

    setToNegativeOneIfNaN(corrShiftedLeft);
    setToNegativeOneIfNaN(corrShiftedRight);

    bool res = corr > 0.7 && corrShiftedLeft - corrWithPadding < 0.1 && corrShiftedRight - corrWithPadding < 0.1;
    return res;
}
}  // namespace CruxQuant
