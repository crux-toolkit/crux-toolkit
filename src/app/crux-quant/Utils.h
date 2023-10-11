#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "IndexedMassSpectralPeak.h"
#include "PpmTolerance.h"
#include "io/MatchFileReader.h"
#include "io/SpectrumCollectionFactory.h"

using std::map;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

namespace CruxQuant {

const int BINS_PER_DALTON = 100;
const double PROTONMASS = 1.007276466879;

const int NUMISOTOPES_REQUIRED = 2;              // May need to make this a user input
const double PEAK_FINDING_PPM_TOLERANCE = 20.0;  // May need to make this a user input
const double PPM_TOLERANCE = 10.0;               // May need to make this a user input
const bool ID_SPECIFIC_CHARGE_STATE = false;  // May need to make this a user input

string calcFormula(string seq);

struct Identification {
    string sequence;
    int charge;
    double peptideMass = 0.0;
    double monoisotopicMass;
    double peakfindingMass;
    double precursorCharge;
    string spectralFile;
    FLOAT_T ms2RetentionTimeInMinutes;
};

struct Ms1ScanInfo {
    int oneBasedScanNumber;
    int zeroBasedMs1ScanIndex;
    double retentionTime;
};

struct IsotopicEnvelope {
    IndexedMassSpectralPeak indexedPeak;
    int chargeState;
    double intensity;

    IsotopicEnvelope() : indexedPeak(),
                         chargeState(0),
                         intensity(0.0) {}
    // Constructor
    IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity)
        : indexedPeak(monoisotopicPeak), chargeState(chargeState), intensity(intensity / chargeState) {}

    void Normalize(double normalizationFactor) {
        intensity *= normalizationFactor;
    }

    bool operator==(const IsotopicEnvelope& otherEnv) const {
        return chargeState == otherEnv.chargeState && indexedPeak == otherEnv.indexedPeak;
    }

    bool operator!=(const IsotopicEnvelope& otherEnv) const {
        return !(*this == otherEnv);
    }

    size_t Hash() const {
        return std::hash<int>()(chargeState) ^ indexedPeak.GetHashCode();
    }

    IsotopicEnvelope& operator=(const IsotopicEnvelope& other) {
    if (this != &other) {
        IsotopicEnvelope temp(other); // Create a copy using the copy constructor
        std::swap(*this, temp);       // Swap the contents with the temporary object
    }
    return *this;
}

};

struct ChromatographicPeak {
    vector<Identification> identifications;
    string spectralFile;
    vector<IsotopicEnvelope> isotopicEnvelopes;
    double intensity;
    double massError;
    int numChargeStatesObserved;
    IsotopicEnvelope apex;
    bool isMbrPeak;
    int _index = std::numeric_limits<int>::max();

    ChromatographicPeak(const Identification& id, bool _isMbrPeak, const string& _spectralFile)
        : isMbrPeak(_isMbrPeak), spectralFile(_spectralFile) {
        isotopicEnvelopes.push_back(IsotopicEnvelope());
        identifications.push_back(id);
        massError = std::numeric_limits<double>::quiet_NaN();
        numChargeStatesObserved = 0;
        _index = std::numeric_limits<int>::max();
    }
    void calculateIntensityForThisFeature(bool integrate);
};

struct IndexedSpectralResults {
    map<int, map<int, IndexedMassSpectralPeak>> _indexedPeaks;

    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;
};

Crux::SpectrumCollection* loadSpectra(const string& file, int ms_level);

IndexedSpectralResults indexedMassSpectralPeaks(Crux::SpectrumCollection* spectrum_collection, const string& spectra_file);

vector<Identification> createIdentifications(MatchFileReader* matchFileReader, const string& spectra_file, Crux::SpectrumCollection* spectrum_collection);

unordered_map<string, vector<pair<double, double>>> calculateTheoreticalIsotopeDistributions(const vector<Identification>& allIdentifications);

void setPeakFindingMass(
    vector<Identification>& allIdentifications, 
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution
);

vector<double> createChargeStates(const vector<Identification>& allIdentifications);

void quantifyMs2IdentifiedPeptides(
    string spectraFile, 
    const vector<Identification>& allIdentifications, 
    const vector<double>& chargeStates
);

double toMz(double mass, int charge);

double toMass(double massToChargeRatio, int charge);

IndexedMassSpectralPeak* getIndexedPeak(
    double theorMass, 
    int zeroBasedScanIndex, 
    PpmTolerance tolerance, 
    int chargeState, 
    map<int, map<int, IndexedMassSpectralPeak>> indexedPeaks
);

void processRange(
    int start, 
    int end, 
    const vector<Identification>& ms2IdsForThisFile, 
    const string& spectralFile,
    const vector<double>& chargeStates,
    vector<ChromatographicPeak>& chromatographicPeaks,
    PpmTolerance& peakfindingTol
);

}  // namespace CruxQuant
