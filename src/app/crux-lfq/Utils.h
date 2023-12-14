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
#include "ProteinGroup.h"
#include "io/carp.h"
#include "util/Params.h"


using std::map;
using std::pair;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

typedef std::tuple<double, double, double> IsotopePeak;

namespace CruxQuant {

const int BINS_PER_DALTON = 100;
const double PROTONMASS = 1.007276466879;
const double C13MinusC12 = 1.00335483810;

extern int NUM_ISOTOPES_REQUIRED;                    // Default value is 2
extern double PEAK_FINDING_PPM_TOLERANCE;         // Default value is 20.0
extern double PPM_TOLERANCE;                      // Default value is 10.0
extern bool ID_SPECIFIC_CHARGE_STATE;            // Default value is false
extern int  MISSED_SCANS_ALLOWED;                    // Default value is 1
extern double ISOTOPE_TOLERANCE_PPM;               // Default value is 5.0
extern bool INTEGRATE;                           // Default value is false
extern double DISCRIMINATION_FACTOR_TO_CUT_PEAK;   // Default value is 0.6
extern bool QUANTIFY_AMBIGUOUS_PEPTIDES;         // Default value is false
extern bool USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT;    // Default value is false
extern bool NORMALIZE;                          // Default value is false
// MBR settings
extern bool MATCH_BETWEEN_RUNS; // Default value is false
extern double MATCH_BETWEEN_RUNS_PPM_TOLERANCE; // Default value is 10.0
extern double MAX_MBR_WINDOW; // Default value is 2.5
extern bool REQUIRE_MSMS_ID_IN_CONDITION; // Default value is false

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
    bool useForProteinQuant;
    unordered_set<ProteinGroup> proteinGroups;


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
            && modifications == other.modifications;
            // && posteriorErrorProbability == other.posteriorErrorProbability;
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

class Peptide;

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

enum class DetectionType {
        MSMS,
        MBR,
        NotDetected,
        MSMSAmbiguousPeakfinding,
        MSMSIdentifiedButNotQuantified,
        Imputed
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

namespace std {
    template <>
    struct hash<CruxQuant::Identification>{
        size_t operator()(const CruxQuant::Identification& id) const {
            return 
                hash<string>()(id.sequence) ^ 
                hash<int>()(id.charge) ^ 
                hash<double>()(id.peptideMass) ^
                hash<double>()(id.monoIsotopicMass) ^
                hash<double>()(id.peakFindingMass) ^
                hash<double>()(id.precursorCharge) ^
                hash<string>()(id.spectralFile) ^
                hash<double>()(id.ms2RetentionTimeInMinutes) ^
                hash<int>()(id.scanId) ^
                hash<string>()(id.modifications);
        }
    };

    template<>
    struct equal_to<CruxQuant::Identification> {
        bool operator()(const CruxQuant::Identification& id1, const CruxQuant::Identification& id2) const {
            // Compare id1 and id2 for equality
            // This uses the operator== that is already defined
            return id1 == id2;
        }
    };
}