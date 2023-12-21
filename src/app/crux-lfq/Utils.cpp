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
#include <tuple>
#include <unordered_map>
#include <vector>

#include "CMercury8.h"
#include "CQStatistics.h"
#include "ChromatographicPeak.h"
#include "Results.h"
#include "csv.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

std::mutex mtx;  // Declare a mutex

const int MaxThreads = 4;

namespace CruxLFQ {

map<int, PSM> create_psm_map(const string& psm_file,
                             const string& psm_file_format) {
    carp(CARP_INFO, "loading psm data ...");
    map<int, PSM> psm_datum;

    string sequence_col, modifications_col;
    int scan_col, charge_col;
    double peptide_mass_col, spectrum_precursor_mz_col;

    // Even though the headers for assign-confidence and percolator are the same,
    // the code is split for the purposes of seperation of concerns

    if (psm_file_format == "tide-search") {
        io::CSVReader<6, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>>
            matchFileReader(psm_file);
        matchFileReader.read_header(io::ignore_extra_column, "scan", "charge",
                                    "spectrum precursor m/z", "peptide mass",
                                    "sequence", "modifications");
        while (matchFileReader.read_row(scan_col, charge_col,
                                        spectrum_precursor_mz_col, peptide_mass_col,
                                        sequence_col, modifications_col)) {
            PSM psm = {sequence_col,
                       scan_col,
                       charge_col,
                       peptide_mass_col,
                       spectrum_precursor_mz_col,
                       modifications_col};
            psm_datum[scan_col] = psm;
        }
    } else if (psm_file_format == "assign-confidence") {
        io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>>
            matchFileReader(psm_file);
        matchFileReader.read_header(io::ignore_extra_column, "scan", "charge",
                                    "spectrum precursor m/z", "peptide mass",
                                    "sequence");
        string default_modifications_col = "";
        while (matchFileReader.read_row(scan_col, charge_col,
                                        spectrum_precursor_mz_col, peptide_mass_col,
                                        sequence_col)) {
            PSM psm = {sequence_col,
                       scan_col,
                       charge_col,
                       peptide_mass_col,
                       spectrum_precursor_mz_col,
                       default_modifications_col};
            psm_datum[scan_col] = psm;
        }
    } else if (psm_file_format == "percolator") {
        io::CSVReader<5, io::trim_chars<' ', '\t'>, io::no_quote_escape<'\t'>>
            matchFileReader(psm_file);
        matchFileReader.read_header(io::ignore_extra_column, "scan", "charge",
                                    "spectrum precursor m/z", "peptide mass",
                                    "sequence");
        string default_modifications_col = "";
        while (matchFileReader.read_row(scan_col, charge_col,
                                        spectrum_precursor_mz_col, peptide_mass_col,
                                        sequence_col)) {
            PSM psm = {sequence_col,
                       scan_col,
                       charge_col,
                       peptide_mass_col,
                       spectrum_precursor_mz_col,
                       default_modifications_col};
            psm_datum[scan_col] = psm;
        }

        return psm_datum;
    }
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
                if (x == 147) O += 1;
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

unordered_map<string, vector<pair<double, double>>>
calculateTheoreticalIsotopeDistributions(
    const vector<Identification>& allIdentifications) {
    unordered_map<string, vector<pair<double, double>>>
        modifiedSequenceToIsotopicDistribution;

    for (const auto& identification : allIdentifications) {
        string peptide_sequence = identification.sequence;
        int charge = identification.precursorCharge;
        double peptide_mass = identification.peptideMass;

        if (modifiedSequenceToIsotopicDistribution.find(peptide_sequence) !=
            modifiedSequenceToIsotopicDistribution.end()) {
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

        double highestAbundance =
            *std::max_element(abundances.begin(), abundances.end());

        for (int i = 0; i < masses.size(); i++) {
            // Calculate the expected isotopic mass shift for this peptide
            masses[i] -= peptide_mass;

            // Normalize the abundance of each isotope
            abundances[i] /= highestAbundance;

            // Look for these isotopes
            if (isotopicMassesAndNormalizedAbundances.size() <
                    NUM_ISOTOPES_REQUIRED ||
                abundances[i] > 0.1) {
                isotopicMassesAndNormalizedAbundances.push_back(
                    std::make_pair(masses[i], abundances[i]));
            }
        }
        modifiedSequenceToIsotopicDistribution[peptide_sequence] =
            isotopicMassesAndNormalizedAbundances;
        delete[] char_array;
    }

    std::unordered_map<std::string, std::vector<Identification>>
        peptideModifiedSequences;

    for (const auto& identification : allIdentifications) {
        peptideModifiedSequences[identification.sequence].push_back(identification);
    }

    for (auto& identifications : peptideModifiedSequences) {
        // isotope where normalized abundance is 1
        double mostAbundantIsotopeShift = 0.0;

        auto& isotopicDistribution =
            modifiedSequenceToIsotopicDistribution[identifications.first];

        for (const auto& item : isotopicDistribution) {
            if (item.second == 1.0) {
                mostAbundantIsotopeShift = item.first;
                break;
            }
        }

        for (auto& identification : identifications.second) {
            identification.peakFindingMass =
                identification.monoIsotopicMass + mostAbundantIsotopeShift;
        }
    }

    return modifiedSequenceToIsotopicDistribution;
}

void setPeakFindingMass(vector<Identification>& allIdentifications,
                        unordered_map<string, vector<pair<double, double>>>&
                            modifiedSequenceToIsotopicDistribution) {
    for (auto& identification : allIdentifications) {
        const string& sequence = identification.sequence;

        // Find the isotope where normalized abundance is 1
        double mostAbundantIsotopeShift = 0.0;

        const auto& isotopicDistribution =
            modifiedSequenceToIsotopicDistribution.find(sequence);
        if (isotopicDistribution != modifiedSequenceToIsotopicDistribution.end()) {
            for (const auto& item : isotopicDistribution->second) {
                if (item.second == 1.0) {
                    mostAbundantIsotopeShift = item.first;
                    break;
                }
            }
        }

        identification.peakFindingMass =
            identification.monoIsotopicMass + mostAbundantIsotopeShift;
    }
}

vector<int> createChargeStates(
    const vector<Identification>& allIdentifications) {
    auto minChargeState = std::numeric_limits<int>::max();
    auto maxChargeState = std::numeric_limits<int>::lowest();
    for (const auto& identification : allIdentifications) {
        minChargeState = std::min(minChargeState, identification.precursorCharge);
        maxChargeState = std::max(maxChargeState, identification.precursorCharge);
    }
    vector<int> chargeStates;

    for (int charge = minChargeState; charge <= maxChargeState; charge++) {
        chargeStates.push_back(charge);
    }

    return chargeStates;
}

void processRange(int start, int end,
                  const vector<Identification>& ms2IdsForThisFile,
                  const string& spectralFile, const vector<int>& chargeStates,
                  vector<ChromatographicPeak>& chromatographicPeaks,
                  PpmTolerance& peakfindingTol,
                  unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
                  map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
                  PpmTolerance& ppmTolerance,
                  unordered_map<string, vector<pair<double, double>>>&
                      modifiedSequenceToIsotopicDistribution,
                  CruxLFQResults& lfqResults) {
    std::lock_guard<std::mutex> lock(mtx);
    for (int i = start; i < end; ++i) {
        const Identification& identification = ms2IdsForThisFile[i];

        ChromatographicPeak msmsFeature(identification, false, spectralFile);

        // msmsFeature._index = i;
        chromatographicPeaks.push_back(msmsFeature);
        for (const auto& chargeState : chargeStates) {
            if (ID_SPECIFIC_CHARGE_STATE &&
                chargeState != identification.precursorCharge) {
                continue;
            }

            vector<IndexedMassSpectralPeak*> xic = peakFind(
                identification.ms2RetentionTimeInMinutes,  // This value may not be in minutes, check and convert
                identification.peakFindingMass, chargeState,
                identification.spectralFile, peakfindingTol, _ms1Scans,
                indexedPeaks);

            // Sort the xic vector based on the RetentionTime
            std::sort(xic.begin(), xic.end(),
                      [](const IndexedMassSpectralPeak* a,
                         const IndexedMassSpectralPeak* b) {
                          return a->retentionTime < b->retentionTime;
                      });

            xic.erase(std::remove_if(xic.begin(), xic.end(),
                                     [&](const IndexedMassSpectralPeak* p) {
                                         return !ppmTolerance.Within(
                                             toMass(p->mz, chargeState),
                                             identification.peakFindingMass);
                                     }),
                      xic.end());

            vector<IsotopicEnvelope> isotopicEnvelopes = getIsotopicEnvelopes(
                xic, identification, chargeState,
                modifiedSequenceToIsotopicDistribution, indexedPeaks);

            msmsFeature.isotopicEnvelopes.insert(
                msmsFeature.isotopicEnvelopes.end(), isotopicEnvelopes.begin(),
                isotopicEnvelopes.end());
        }

        msmsFeature.calculateIntensityForThisFeature(INTEGRATE);

        cutPeak(msmsFeature, identification.ms2RetentionTimeInMinutes, _ms1Scans);

        if (msmsFeature.isotopicEnvelopes.empty()) {
            continue;
        }

        vector<IsotopicEnvelope> precursorXic;
        std::copy_if(msmsFeature.isotopicEnvelopes.begin(),
                     msmsFeature.isotopicEnvelopes.end(),
                     std::back_inserter(precursorXic),
                     [identification](const IsotopicEnvelope& p) {
                         return p.chargeState == identification.precursorCharge;
                     }

        );

        if (precursorXic.empty()) {
            msmsFeature.isotopicEnvelopes.clear();
            continue;
        }
        // Find the minimum and maximum values
        int min =
            std::numeric_limits<int>::max();  // Initialize with a large value
        int max =
            std::numeric_limits<int>::min();  // Initialize with a small value

        for (const IsotopicEnvelope& p : precursorXic) {
            int scanIndex = p.indexedPeak.zeroBasedMs1ScanIndex;
            if (scanIndex < min) {
                min = scanIndex;
            }
            if (scanIndex > max) {
                max = scanIndex;
            }
        }

        msmsFeature.isotopicEnvelopes.erase(
            std::remove_if(msmsFeature.isotopicEnvelopes.begin(),
                           msmsFeature.isotopicEnvelopes.end(),
                           [min](const IsotopicEnvelope& p) {
                               return p.indexedPeak.zeroBasedMs1ScanIndex < min;
                           }),
            msmsFeature.isotopicEnvelopes.end());

        msmsFeature.isotopicEnvelopes.erase(
            std::remove_if(msmsFeature.isotopicEnvelopes.begin(),
                           msmsFeature.isotopicEnvelopes.end(),
                           [max](const IsotopicEnvelope& p) {
                               return p.indexedPeak.zeroBasedMs1ScanIndex > max;
                           }),
            msmsFeature.isotopicEnvelopes.end());

        msmsFeature.calculateIntensityForThisFeature(INTEGRATE);
    }

    lfqResults.Peaks[spectralFile].insert(lfqResults.Peaks[spectralFile].end(),
                                          chromatographicPeaks.begin(),
                                          chromatographicPeaks.end());
}

void quantifyMs2IdentifiedPeptides(
    string spectraFile,
    const vector<Identification>& allIdentifications,
    const vector<int>& chargeStates,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
    unordered_map<string, vector<pair<double, double>>>&
        modifiedSequenceToIsotopicDistribution,
    CruxLFQResults& lfqResults) {
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
        int end = (i == MaxThreads - 1) ? (start + chunkSize + remainder)
                                        : (start + chunkSize);

        threads.push_back(
            std::thread([start, end, &ms2IdsForThisFile, &spectraFile,
                         &chargeStates, &chromatographicPeaks, &peakfindingTol,
                         &_ms1Scans, &indexedPeaks, &ppmTolerance,
                         &modifiedSequenceToIsotopicDistribution, &lfqResults]() {
                processRange(start, end, ms2IdsForThisFile, spectraFile, chargeStates,
                             chromatographicPeaks, peakfindingTol, _ms1Scans,
                             indexedPeaks, ppmTolerance,
                             modifiedSequenceToIsotopicDistribution, lfqResults);
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

IndexedMassSpectralPeak* getIndexedPeak(
    double theorMass, int zeroBasedScanIndex, PpmTolerance tolerance,
    int chargeState,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
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
        if (j < indexedPeaks.size() && indexedPeaks[j].size() > 0) {
            map<int, IndexedMassSpectralPeak> bin = indexedPeaks[j];
            auto startIterator = bin.find(zeroBasedScanIndex);
            for (auto it = startIterator; it != bin.end(); ++it) {
                auto _peak = it->second;
                if (_peak.zeroBasedMs1ScanIndex > zeroBasedScanIndex) {
                    break;
                }

                double expMass = toMass(_peak.mz, chargeState);

                if (tolerance.Within(expMass, theorMass) &&
                    _peak.zeroBasedMs1ScanIndex == zeroBasedScanIndex &&
                    (bestPeak == nullptr ||
                     std::abs(expMass - theorMass) <
                         std::abs(toMass(bestPeak->mz, chargeState) - theorMass))) {
                    bestPeak = &_peak;
                }
            }
        }
    }
    return bestPeak;
}

vector<IndexedMassSpectralPeak*> peakFind(
    double idRetentionTime, double mass, int charge, const string& spectra_file,
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

    std::sort(
        xic.begin(), xic.end(),
        [](const IndexedMassSpectralPeak* x, const IndexedMassSpectralPeak* y) {
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
    const Identification& identification, const int chargeState,
    unordered_map<string, vector<pair<double, double>>>&
        modifiedSequenceToIsotopicDistribution,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
    vector<IsotopicEnvelope> isotopicEnvelopes;
    vector<pair<double, double>> isotopeMassShifts =
        modifiedSequenceToIsotopicDistribution[identification.sequence];

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

    int peakfindingMassIndex = static_cast<int>(
        round(identification.peakFindingMass - identification.monoIsotopicMass));

    // For each peak in the XIC, we consider the possibility that there was an
    // off-by-one or missed monoisotopic mass error in peak assignment /
    // deconvolution. The -1 key in this dictionary corresponds to a negative
    // off-by-one error, the +1 key corresponds to a positive off-by-one error,
    // and the 0 key corresponds to accurate assignment/deconvolution.
    map<int, vector<IsotopePeak>> massShiftToIsotopePeaks = {
        {-1, {}}, {0, {}}, {1, {}}};

    vector<int> directions = {-1, 1};

    for (const IndexedMassSpectralPeak* peak : xic) {
        std::fill_n(experimentalIsotopeIntensities.begin(),
                    experimentalIsotopeIntensities.size(), 0);

        for (auto& kvp : massShiftToIsotopePeaks) {
            kvp.second.clear();
        }

        // isotope masses are calculated relative to the observed peak
        double observedMass = toMass(peak->mz, chargeState);
        double observedMassError = observedMass - identification.peakFindingMass;

        for (auto& shift : massShiftToIsotopePeaks) {
            // look for each isotope peak in the data
            // This is done by starting with the first isotope with mass less than the
            // peak finding (most abundant) mass, then working backwards to find every
            // isotope with mass < most abundant mass. Once an expected isotopic peak
            // can not be found, the loop breaks and we begin working our way forward,
            // starting with the peak finding mass and locating every peak with mass >
            // peak finding mass. Once an expected isotopic peak can not be found, it
            // is assumed that we have located every isotope present and the loop
            // breaks.

            for (int direction : directions) {
                int start =
                    (direction == -1) ? peakfindingMassIndex - 1 : peakfindingMassIndex;

                for (int i = start; i < theoreticalIsotopeAbundances.size() && i >= 0;
                     i += direction) {
                    double isotopeMass =
                        identification.monoIsotopicMass + observedMassError +
                        theoreticalIsotopeMassShifts[i] + shift.first * C13MinusC12;
                    double theoreticalIsotopeIntensity =
                        theoreticalIsotopeAbundances[i] * peak->intensity;

                    IndexedMassSpectralPeak* isotopePeak =
                        getIndexedPeak(isotopeMass, peak->zeroBasedMs1ScanIndex,
                                       isotopeTolerance, chargeState, indexedPeaks);

                    if (isotopePeak == nullptr ||
                        isotopePeak->intensity < theoreticalIsotopeIntensity / 4.0 ||
                        isotopePeak->intensity > theoreticalIsotopeIntensity * 4.0) {
                        break;
                    }

                    shift.second.push_back(std::make_tuple(isotopePeak->intensity,
                                                           theoreticalIsotopeIntensity,
                                                           isotopeMass));
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
            if (checkIsotopicEnvelopeCorrelation(massShiftToIsotopePeaks, peak,
                                                 chargeState, isotopeTolerance,
                                                 indexedPeaks)) {
                for (size_t i = 0; i < experimentalIsotopeIntensities.size(); ++i) {
                    if (experimentalIsotopeIntensities[i] == 0) {
                        experimentalIsotopeIntensities[i] =
                            theoreticalIsotopeAbundances[i] *
                            experimentalIsotopeIntensities[peakfindingMassIndex];
                    }
                }
            }

            double sumExperimentalIntensities =
                std::accumulate(experimentalIsotopeIntensities.begin(),
                                experimentalIsotopeIntensities.end(), 0.0);

            IsotopicEnvelope isotopicEnvelope(*peak, chargeState,
                                              sumExperimentalIntensities);

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

    filterResults results{expIntensity, theorIntensity};

    return results;
}

void setToNegativeOneIfNaN(double& value) {
    if (std::isnan(value)) {
        value = -1;
    }
}

bool checkIsotopicEnvelopeCorrelation(
    map<int, vector<IsotopePeak>>& massShiftToIsotopePeaks,
    const IndexedMassSpectralPeak* peak, int chargeState,
    PpmTolerance& isotopeTolerance,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks) {
    filterResults results =
        filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[0]);
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
        IndexedMassSpectralPeak* unexpectedPeak =
            getIndexedPeak(unexpectedMass, peak->zeroBasedMs1ScanIndex,
                           isotopeTolerance, chargeState, indexedPeaks);

        if (unexpectedPeak == nullptr) {
            shift.second.push_back(std::make_tuple(0.0, 0.0, unexpectedMass));
        } else {
            shift.second.push_back(
                std::make_tuple(unexpectedPeak->intensity, 0.0, unexpectedMass));
        }
    }

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[0]);
    double corrWithPadding =
        Pearson(results.expIntensity, results.theorIntensity);

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[-1]);
    double corrShiftedLeft =
        Pearson(results.expIntensity, results.theorIntensity);

    results = filterMassShiftToIsotopePeaks(massShiftToIsotopePeaks[1]);
    double corrShiftedRight =
        Pearson(results.expIntensity, results.theorIntensity);

    setToNegativeOneIfNaN(corrShiftedLeft);
    setToNegativeOneIfNaN(corrShiftedRight);

    bool res = corr > 0.7 && corrShiftedLeft - corrWithPadding < 0.1 &&
               corrShiftedRight - corrWithPadding < 0.1;
    return res;
}

void cutPeak(ChromatographicPeak& peak, double identificationTime,
             unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans) {
    // find out if we need to split this peak by using the discrimination factor
    // this method assumes that the isotope envelopes in a chromatographic peak
    // are already sorted by MS1 scan number
    bool cutThisPeak = false;

    if (peak.isotopicEnvelopes.size() < 5) {
        return;
    }

    vector<IsotopicEnvelope> timePointsForApexZ;

    std::copy_if(peak.isotopicEnvelopes.begin(), peak.isotopicEnvelopes.end(),
                 std::back_inserter(timePointsForApexZ),
                 [&peak](const IsotopicEnvelope& p) {
                     return p.chargeState == peak.apex.chargeState;
                 });

    std::unordered_set<int> scanNumbers;
    std::transform(timePointsForApexZ.begin(), timePointsForApexZ.end(),
                   std::inserter(scanNumbers, scanNumbers.end()),
                   [](const IsotopicEnvelope& p) {
                       return p.indexedPeak.zeroBasedMs1ScanIndex;
                   });

    auto it = std::find(timePointsForApexZ.begin(), timePointsForApexZ.end(),
                        peak.apex);
    int apexIndex = -1;  // Initialize to -1 if the element is not found

    if (it != timePointsForApexZ.end()) {
        // Calculate the index by subtracting the iterator from the beginning of the
        // vector
        apexIndex = std::distance(timePointsForApexZ.begin(), it);
    }

    IsotopicEnvelope* valleyEnvelope = nullptr;
    int directions[] = {1, -1};
    for (int direction : directions) {
        valleyEnvelope = nullptr;
        int indexOfValley = 0;
        int lastIndex = timePointsForApexZ.size() - 1;
        for (int i = apexIndex + direction; i <= lastIndex && i >= 0;
             i += direction) {
            IsotopicEnvelope timepoint = timePointsForApexZ[i];

            if (valleyEnvelope == nullptr ||
                timepoint.intensity < valleyEnvelope->intensity) {
                valleyEnvelope = &timepoint;
                indexOfValley = i;
            }

            double discriminationFactor =
                (timepoint.intensity - valleyEnvelope->intensity) /
                timepoint.intensity;

            if (discriminationFactor > DISCRIMINATION_FACTOR_TO_CUT_PEAK &&
                (indexOfValley + direction < timePointsForApexZ.size() &&
                 indexOfValley + direction >= 0)) {
                IsotopicEnvelope secondValleyTimepoint =
                    timePointsForApexZ[indexOfValley + direction];
                discriminationFactor =
                    (timepoint.intensity - secondValleyTimepoint.intensity) /
                    timepoint.intensity;
                if (discriminationFactor > DISCRIMINATION_FACTOR_TO_CUT_PEAK) {
                    cutThisPeak = true;
                    break;
                }
                int nextMs1ScanNum = -1;
                for (int j = valleyEnvelope->indexedPeak.zeroBasedMs1ScanIndex - 1;
                     j < _ms1Scans[peak.spectralFile].size() && j >= 0;
                     j += direction) {
                    if (_ms1Scans[peak.spectralFile][j].oneBasedScanNumber >= 0 &&
                        _ms1Scans[peak.spectralFile][j].oneBasedScanNumber !=
                            valleyEnvelope->indexedPeak.zeroBasedMs1ScanIndex) {
                        nextMs1ScanNum = j + 1;
                        break;
                    }
                }

                if (scanNumbers.find(nextMs1ScanNum) == scanNumbers.end()) {
                    cutThisPeak = true;
                    break;
                }
            }
        }
        if (cutThisPeak) {
            break;
        }
    }
    if (cutPeak) {
        if (identificationTime > valleyEnvelope->indexedPeak.retentionTime) {
            peak.isotopicEnvelopes.erase(
                std::remove_if(peak.isotopicEnvelopes.begin(),
                               peak.isotopicEnvelopes.end(),
                               [&valleyEnvelope](const IsotopicEnvelope p) {
                                   return p.indexedPeak.retentionTime <=
                                          valleyEnvelope->indexedPeak.retentionTime;
                               }),
                peak.isotopicEnvelopes.end());
        } else {
            peak.isotopicEnvelopes.erase(
                std::remove_if(peak.isotopicEnvelopes.begin(),
                               peak.isotopicEnvelopes.end(),
                               [&valleyEnvelope](const IsotopicEnvelope p) {
                                   return p.indexedPeak.retentionTime >=
                                          valleyEnvelope->indexedPeak.retentionTime;
                               }),
                peak.isotopicEnvelopes.end());
        }

        // recalculate intensity for the peak
        peak.calculateIntensityForThisFeature(INTEGRATE);
        peak.SplitRT = valleyEnvelope->indexedPeak.retentionTime;

        // recursively cut
        cutPeak(peak, identificationTime, _ms1Scans);
    }
}

void runErrorChecking(const string& spectraFile, CruxLFQResults& lfqResults) {
    carp(CARP_INFO, "Checking errors");

    lfqResults.Peaks[spectraFile].erase(
        std::remove_if(lfqResults.Peaks[spectraFile].begin(),
                       lfqResults.Peaks[spectraFile].end(),
                       [](const ChromatographicPeak& p) {
                           return &p == nullptr ||
                                  (p.isMbrPeak && p.isotopicEnvelopes.empty());
                       }),
        lfqResults.Peaks[spectraFile].end());

    // merge duplicate peaks and handle MBR/MSMS peakfinding conflicts
    unordered_map<IndexedMassSpectralPeak, ChromatographicPeak>
        errorCheckedPeaksGroupedByApex;
    vector<ChromatographicPeak> errorCheckedPeaks;

    std::sort(lfqResults.Peaks[spectraFile].begin(),
              lfqResults.Peaks[spectraFile].end(),
              [](const ChromatographicPeak& a, const ChromatographicPeak& b) {
                  return a.isMbrPeak < b.isMbrPeak;
              });

    for (ChromatographicPeak& tryPeak : lfqResults.Peaks[spectraFile]) {
        tryPeak.calculateIntensityForThisFeature(INTEGRATE);
        tryPeak.resolveIdentifications();

        IsotopicEnvelope emptyStruct = {};

        if (memcmp(&tryPeak.apex, &emptyStruct, sizeof(IsotopicEnvelope)) == 0) {
            if (tryPeak.isMbrPeak) {
                continue;
            }
            errorCheckedPeaks.push_back(tryPeak);
            continue;
        }

        IndexedMassSpectralPeak apexImsPeak = tryPeak.apex.indexedPeak;
        auto it = errorCheckedPeaksGroupedByApex.find(apexImsPeak);
        if (it != errorCheckedPeaksGroupedByApex.end()) {
            ChromatographicPeak storedPeak = it->second;
            if (tryPeak.isMbrPeak && &storedPeak == nullptr) {
                continue;
            }
            if (!tryPeak.isMbrPeak && !storedPeak.isMbrPeak) {
                storedPeak.mergeFeatureWith(tryPeak, INTEGRATE);
            } else if (tryPeak.isMbrPeak && !storedPeak.isMbrPeak) {
                continue;
            } else if (tryPeak.isMbrPeak && storedPeak.isMbrPeak) {
                // Check if the ModifiedSequence of the first Identification of tryPeak
                // and storedPeak is the same
                if (tryPeak.identifications.front().modifications ==
                    storedPeak.identifications.front().modifications) {
                    // Merge the features if the ModifiedSequences match
                    storedPeak.mergeFeatureWith(tryPeak, INTEGRATE);
                } else if (tryPeak.MbrScore >
                           storedPeak
                               .MbrScore) {  // The else if block is redundant since in
                                             // our use case MbrScore will not be set
                    // Compare MbrScores and update the errorCheckedPeaksGroupedByApex if
                    // necessary
                    errorCheckedPeaksGroupedByApex.insert(
                        std::make_pair(tryPeak.apex.indexedPeak, tryPeak));
                }
            }

        } else {
            errorCheckedPeaksGroupedByApex.insert({apexImsPeak, tryPeak});
        }
    }
    // Assuming 'errorCheckedPeaks' is a std::vector and
    // 'errorCheckedPeaksGroupedByApex' is a std::map

    for (const auto& pair : errorCheckedPeaksGroupedByApex) {
        const ChromatographicPeak& value = pair.second;
        if (&value != nullptr) {
            errorCheckedPeaks.push_back(value);
        }
    }

    lfqResults.Peaks[spectraFile] = errorCheckedPeaks;
}

// void quantifyMatchBetweenRunsPeaks(const string& spectraFile, CruxLFQResults&
// lfqResults) {
// }
}  // namespace CruxLFQ
