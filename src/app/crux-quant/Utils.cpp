#pragma once

#include "Utils.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "CMercury8.h"
#include "PpmTolerance.h"

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;



namespace CruxQuant {

Crux::SpectrumCollection* loadSpectra(const string &file, int ms_level) {
    Crux::SpectrumCollection *spectra = SpectrumCollectionFactory::create(file);
    spectra->parse(ms_level = ms_level);
    return spectra;
}

IndexedSpectralResults indexedMassSpectralPeaks(Crux::SpectrumCollection *spectrum_collection, const string &spectra_file) {
    carp(CARP_INFO, "Read %d spectra. for MS1", spectrum_collection->getNumSpectra());

    string _spectra_file(spectra_file);
    // Define the hashmap with a vector of IndexedMassSpectralPeak instances as the value type
    unordered_map<int, vector<IndexedMassSpectralPeak>> _indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;
    _ms1Scans[_spectra_file] = vector<Ms1ScanInfo>();

    IndexedSpectralResults index_results{_indexedPeaks,_ms1Scans};

    if (spectrum_collection->getNumSpectra() <= 0) {
        return index_results;
    }

    int scanIndex = 0;
    int oneBasedScanNumber = 1;

    for (auto spectrum = spectrum_collection->begin(); spectrum != spectrum_collection->end(); ++spectrum) {
        if (*spectrum != nullptr) {
            for (auto peak = (*spectrum)->begin(); peak != (*spectrum)->end(); ++peak) {
                if (*peak != nullptr) {
                    FLOAT_T mz = (*peak)->getLocation();
                    int roundedMz = static_cast<int>(std::round(mz * BINS_PER_DALTON));
                    double retentionTime = (*spectrum)->getRetentionTime();  // retentionTime
                    IndexedMassSpectralPeak spec_data(
                        mz,                              // mz value
                        (*peak)->getIntensity(),         // intensity
                        scanIndex,                       // zeroBasedMs1ScanIndex
                        retentionTime
                    );

                    // Find the corresponding entry in _indexedPeaks and create a new entry if it doesn't exist
                    auto it = _indexedPeaks.find(roundedMz);
                    if (it == _indexedPeaks.end()) {
                        _indexedPeaks.insert({roundedMz, vector<IndexedMassSpectralPeak>()});
                    } else {
                        // Add a new IndexedMassSpectralPeak object to the vector at the corresponding roundedMz entry
                        _indexedPeaks[roundedMz].emplace_back(spec_data);
                    }
                    Ms1ScanInfo scan = {oneBasedScanNumber, scanIndex, retentionTime};
                    _ms1Scans[spectra_file].push_back(scan);
                }
            }
            scanIndex++;
            oneBasedScanNumber++;
        }
    }

    return index_results; 
}

vector<Identification> createIdentifications(MatchFileReader* matchFileReader, const string &spectra_file) {
    vector<Identification> allIdentifications;
    string _spectra_file(spectra_file);
    while (matchFileReader->hasNext()) {
        Identification identification;

        identification.Sequence = matchFileReader->getString(SEQUENCE_COL);
        identification.MonoisotopicMass = matchFileReader->getDouble(PEPTIDE_MASS_COL);
        identification.Charge = matchFileReader->getInteger(CHARGE_COL);
        identification.PeptideMass = matchFileReader->getDouble(PEPTIDE_MASS_COL);
        identification.PrecursorCharge = matchFileReader->getDouble(SPECTRUM_PRECURSOR_MZ_COL);
        identification.spectralFile = _spectra_file;
        allIdentifications.push_back(identification);

        matchFileReader->next();
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
        string peptide_sequence = identification.Sequence;
        int charge = identification.Charge;
        double peptide_mass = identification.PeptideMass;

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
            if (isotopicMassesAndNormalizedAbundances.size() < NUMISOTOPES_REQUIRED || abundances[i] > 0.1) {
                isotopicMassesAndNormalizedAbundances.push_back(std::make_pair(masses[i], abundances[i]));
            }
        }
        modifiedSequenceToIsotopicDistribution[peptide_sequence] = isotopicMassesAndNormalizedAbundances;
        delete[] char_array;
    }

    return modifiedSequenceToIsotopicDistribution;
}

void SetPeakFindingMass(vector<Identification>& allIdentifications, unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution) {
    for (auto& identification : allIdentifications) {
        const string& sequence = identification.Sequence;

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

        identification.PeakfindingMass = identification.MonoisotopicMass + mostAbundantIsotopeShift;
    }
}

vector<double> createChargeStates(const vector<Identification>& allIdentifications) {
    auto minChargeState = std::numeric_limits<double>::max();
    auto maxChargeState = std::numeric_limits<double>::lowest();
    for (const auto& identification : allIdentifications) {
        minChargeState = std::min(minChargeState, identification.PrecursorCharge);
        maxChargeState = std::max(maxChargeState, identification.PrecursorCharge);
    }
    vector<double> chargeStates;

    for (double charge = minChargeState; charge <= maxChargeState; charge++) {
        chargeStates.push_back(charge);
    }

    return chargeStates;
}

void QuantifyMs2IdentifiedPeptides(string spectraFile, const vector<Identification>& allIdentifications){
    vector<Identification> ms2IdsForThisFile;

    // Use std::copy_if to filter the identifications
    std::copy_if(
        allIdentifications.begin(),
        allIdentifications.end(),
        std::back_inserter(ms2IdsForThisFile),
        [&spectraFile](const Identification& id) {
            return id.spectralFile == spectraFile;
        }
    );

    // Check if ms2IdsForThisFile is empty
    if (ms2IdsForThisFile.empty()) {
        return; // Early return
    }

    PpmTolerance peakfindingTol(CruxQuant::PEAK_FINDING_PPM_TOLERANCE); // Peak finding tolerance is generally higher than ppmTolerance
    PpmTolerance ppmTolerance(CruxQuant::PPM_TOLERANCE);

}


}  // namespace CruxQuant
