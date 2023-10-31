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
#include "io/carp.h"


using std::map;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

typedef std::tuple<double, double, double> IsotopePeak;

namespace CruxQuant {

const int BINS_PER_DALTON = 100;
const double PROTONMASS = 1.007276466879;
const double C13MinusC12 = 1.00335483810;

const int NUM_ISOTOPES_REQUIRED = 2;                    // May need to make this a user input
const double PEAK_FINDING_PPM_TOLERANCE = 20.0;         // May need to make this a user input
const double PPM_TOLERANCE = 10.0;                      // May need to make this a user input
const bool ID_SPECIFIC_CHARGE_STATE = false;            // May need to make this a user input
const int  MISSED_SCANS_ALLOWED = 1;                    // May need to make this a user input
const double ISOTOPE_TOLERANCE_PPM = 5.0;               // May need to make this a user input
const bool INTEGRATE = false;                           // May need to make this a user input
const double DISCRIMINATION_FACTOR_TO_CUT_PEAK = 0.6;   // May need to make this a user input

string calcFormula(string seq);

struct Identification {
    string sequence;
    int charge;
    double peptideMass = 0.0;
    double monoIsotopicMass;
    double peakFindingMass;
    double precursorCharge;
    string spectralFile;
    FLOAT_T ms2RetentionTimeInMinutes;
    int scanId;
    string modifications;
    double posteriorErrorProbability = 0; //This may be removed cos it's redundant

    bool operator==(const Identification& other) const {
        // Implement the logic for comparison
        return sequence == other.sequence
            && charge == other.charge
            && peptideMass == other.peptideMass
            && monoIsotopicMass == other.monoIsotopicMass
            && peakFindingMass == other.peakFindingMass
            && precursorCharge == other.precursorCharge
            && spectralFile == other.spectralFile
            && ms2RetentionTimeInMinutes == other.ms2RetentionTimeInMinutes
            && scanId == other.scanId
            && modifications == other.modifications
            && posteriorErrorProbability == other.posteriorErrorProbability;
    }
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
        indexedPeak = other.indexedPeak;
        chargeState = other.chargeState;
        intensity = other.intensity;
    }
    return *this;
}

};

struct ChromatographicPeak;

struct IndexedSpectralResults {
    map<int, map<int, IndexedMassSpectralPeak>> _indexedPeaks;

    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;
};

struct PSM{
    string sequence_col;
    int scan_col;
    int charge_col;
    double peptide_mass_col;
    double spectrum_precursor_mz_col;
    string modifications;
};

unordered_map<string, vector<pair<double, double>>> calculateTheoreticalIsotopeDistributions(const vector<Identification>& allIdentifications);

void setPeakFindingMass(
    vector<Identification>& allIdentifications, 
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution
);

vector<double> createChargeStates(const vector<Identification>& allIdentifications);

// Forward declaration for a class from results.h
class CruxLFQResults; 

void quantifyMs2IdentifiedPeptides(
    string spectraFile, 
    const vector<Identification>& allIdentifications, 
    const vector<double>& chargeStates,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution,
    CruxLFQResults& lfqResults
);

double toMz(double mass, int charge);

double toMass(double massToChargeRatio, int charge);

IndexedMassSpectralPeak* getIndexedPeak(
    double theorMass, 
    int zeroBasedScanIndex, 
    PpmTolerance tolerance, 
    int chargeState, 
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks
);

void processRange(
    int start, 
    int end, 
    const vector<Identification>& ms2IdsForThisFile, 
    const string& spectralFile,
    const vector<double>& chargeStates,
    vector<ChromatographicPeak>& chromatographicPeaks,
    PpmTolerance& peakfindingTol,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
    PpmTolerance& ppmTolerance,
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution,
    CruxLFQResults& lfqResults
);

vector<IndexedMassSpectralPeak*> peakFind(
        double idRetentionTime, 
        double mass, 
        int charge, 
        const string& spectra_file, 
        PpmTolerance tolerance,
        unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
        map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks
);

string getScanID(string spectrum_id);

vector<IsotopicEnvelope> getIsotopicEnvelopes(
    const vector<IndexedMassSpectralPeak*>& xic, 
    const Identification& identification,
    const int chargeState, 
    unordered_map<string, vector<pair<double, double>>>& modifiedSequenceToIsotopicDistribution,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks
);

bool checkIsotopicEnvelopeCorrelation(
    map<int, vector<IsotopePeak>>& massShiftToIsotopePeaks,
    const IndexedMassSpectralPeak* peak,
    int chargeState,
    PpmTolerance& isotopeTolerance,
    map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks
);

struct filterResults{
    vector<double> expIntensity;
    vector<double> theorIntensity;
};

filterResults filterMassShiftToIsotopePeaks(vector<IsotopePeak>& isotopePeaks);

void setToNegativeOneIfNaN(double& value);

map<int, PSM> create_psm_map(const string& psm_file);

void cutPeak(ChromatographicPeak& peak, double identificationTime, unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans);

void runErrorChecking(const string& spectraFile, CruxLFQResults& lfqResults);

}  // namespace CruxQuant
