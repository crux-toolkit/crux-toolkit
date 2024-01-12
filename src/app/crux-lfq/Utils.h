#pragma once
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <ostream>
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

namespace CruxLFQ {

const int BINS_PER_DALTON = 100;
const double PROTONMASS = 1.007276466879;
const double C13MinusC12 = 1.00335483810;

extern int NUM_ISOTOPES_REQUIRED;                   // Default value is 2
extern double PEAK_FINDING_PPM_TOLERANCE;           // Default value is 20.0
extern double PPM_TOLERANCE;                        // Default value is 10.0
extern bool ID_SPECIFIC_CHARGE_STATE;               // Default value is false
extern int MISSED_SCANS_ALLOWED;                    // Default value is 1
extern double ISOTOPE_TOLERANCE_PPM;                // Default value is 5.0
extern bool INTEGRATE;                              // Default value is false
extern double DISCRIMINATION_FACTOR_TO_CUT_PEAK;    // Default value is 0.6
extern bool QUANTIFY_AMBIGUOUS_PEPTIDES;            // Default value is false
extern bool USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT;  // Default value is false
extern bool NORMALIZE;                              // Default value is false

string calcFormula(string seq);

struct Identification {
    string sequence;
    double peptideMass;
    double monoIsotopicMass;
    double peakFindingMass = 0.0;
    int precursorCharge;
    string spectralFile;
    double ms2RetentionTimeInMinutes;
    int scanId;
    string modifications;
    double posteriorErrorProbability =
        0;  // This may be removed cos it's redundant
    bool useForProteinQuant = true;
    unordered_set<ProteinGroup> proteinGroups;

    Identification()
        : peptideMass(0), monoIsotopicMass(0), peakFindingMass(0), precursorCharge(0), ms2RetentionTimeInMinutes(0), scanId(0), posteriorErrorProbability(0), useForProteinQuant(false) {}

    bool operator==(const Identification& other) const {
        // Implement the logic for comparison
        return sequence == other.sequence && peptideMass == other.peptideMass &&
               monoIsotopicMass == other.monoIsotopicMass &&
               peakFindingMass == other.peakFindingMass &&
               precursorCharge == other.precursorCharge &&
               spectralFile == other.spectralFile &&
               ms2RetentionTimeInMinutes == other.ms2RetentionTimeInMinutes &&
               scanId == other.scanId && modifications == other.modifications && posteriorErrorProbability == other.posteriorErrorProbability &&
               useForProteinQuant == other.useForProteinQuant &&
               proteinGroups == other.proteinGroups;
    }

    Identification& operator=(const Identification& other) {
        if (this != &other) {
            sequence = other.sequence;
            peptideMass = other.peptideMass;
            monoIsotopicMass = other.monoIsotopicMass;
            peakFindingMass = other.peakFindingMass;
            precursorCharge = other.precursorCharge;
            spectralFile = other.spectralFile;
            ms2RetentionTimeInMinutes = other.ms2RetentionTimeInMinutes;
            scanId = other.scanId;
            modifications = other.modifications;
            posteriorErrorProbability = other.posteriorErrorProbability;
            useForProteinQuant = other.useForProteinQuant;
            proteinGroups = other.proteinGroups;
        }
        return *this;
    }

    Identification(const Identification& other) {
        sequence = other.sequence;
        peptideMass = other.peptideMass;
        monoIsotopicMass = other.monoIsotopicMass;
        peakFindingMass = other.peakFindingMass;
        precursorCharge = other.precursorCharge;
        spectralFile = other.spectralFile;
        ms2RetentionTimeInMinutes = other.ms2RetentionTimeInMinutes;
        scanId = other.scanId;
        modifications = other.modifications;
        posteriorErrorProbability = other.posteriorErrorProbability;
        useForProteinQuant = other.useForProteinQuant;
        proteinGroups = other.proteinGroups;
    }

    Identification(Identification&& other) noexcept
        : sequence(std::move(other.sequence)),
          peptideMass(other.peptideMass),
          monoIsotopicMass(other.monoIsotopicMass),
          peakFindingMass(other.peakFindingMass),
          precursorCharge(other.precursorCharge),
          spectralFile(std::move(other.spectralFile)),
          ms2RetentionTimeInMinutes(other.ms2RetentionTimeInMinutes),
          scanId(other.scanId),
          modifications(std::move(other.modifications)),
          posteriorErrorProbability(other.posteriorErrorProbability),
          useForProteinQuant(other.useForProteinQuant),
          proteinGroups(std::move(other.proteinGroups)) {
        // Reset other's members to default values, if necessary
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

    IsotopicEnvelope() : indexedPeak(), chargeState(0), intensity(0.0) {}
    // Constructor
    IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState,
                     double intensity)
        : indexedPeak(monoisotopicPeak),
          chargeState(chargeState),
          intensity(intensity / chargeState) {}

    void Normalize(double normalizationFactor) {
        intensity *= normalizationFactor;
    }

    bool operator==(const IsotopicEnvelope& otherEnv) const {
        return chargeState == otherEnv.chargeState &&
               indexedPeak == otherEnv.indexedPeak;
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

struct PSM {
    string sequence_col;
    int scan_col;
    int charge_col;
    double peptide_mass_col;
    double monoisotopic_mass_col;
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

unordered_map<string, vector<pair<double, double>>>
calculateTheoreticalIsotopeDistributions(vector<Identification>& allIdentifications);

// void setPeakFindingMass(vector<Identification>& allIdentifications,
// unordered_map<string, vector<pair<double, double>>>&
// modifiedSequenceToIsotopicDistribution);

vector<int> createChargeStates(
    const vector<Identification>& allIdentifications);

// Forward declaration for a class from results.h
class CruxLFQResults;

// void quantifyMs2IdentifiedPeptides(
//     string spectraFile, const vector<Identification>& allIdentifications,
//     const vector<int>& chargeStates,
//     unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
//     const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
//     unordered_map<string, vector<pair<double, double>>>&
//         modifiedSequenceToIsotopicDistribution,
//     CruxLFQResults& lfqResults);

void quantifyMs2IdentifiedPeptides(
    string spectraFile,
    const vector<Identification>& allIdentifications,
    const vector<int>& chargeStates,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
    unordered_map<string, vector<pair<double, double>>>&
        modifiedSequenceToIsotopicDistribution,
    CruxLFQResults& lfqResults);

double toMz(double mass, int charge);

double toMass(double massToChargeRatio, int charge);

IndexedMassSpectralPeak* getIndexedPeak(
    const double& theorMass, int zeroBasedScanIndex, PpmTolerance tolerance,
    int chargeState, map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks);

// void processRange(int start, int end,
//                   const vector<Identification>& ms2IdsForThisFile,
//                   const string& spectralFile, const vector<int>& chargeStates,
//                   vector<ChromatographicPeak>& chromatographicPeaks,
//                   PpmTolerance& peakfindingTol,
//                   unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
//                   const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
//                   PpmTolerance& ppmTolerance,
//                   unordered_map<string, vector<pair<double, double>>>&
//                       modifiedSequenceToIsotopicDistribution,
//                   CruxLFQResults& lfqResults);

void processRange(int start, int end,
                  const vector<Identification>& ms2IdsForThisFile,
                  const string& spectralFile, const vector<int>& chargeStates,
                  vector<ChromatographicPeak>& chromatographicPeaks,
                  PpmTolerance& peakfindingTol,
                  unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
                  const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks,
                  PpmTolerance& ppmTolerance,
                  unordered_map<string, vector<pair<double, double>>>&
                      modifiedSequenceToIsotopicDistribution,
                  CruxLFQResults& lfqResults);

vector<IndexedMassSpectralPeak*> peakFind(
    double idRetentionTime, const double& mass, int charge, const string& spectra_file,
    PpmTolerance tolerance,
    unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans,
    const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks);

string getScanID(string spectrum_id);

vector<IsotopicEnvelope> getIsotopicEnvelopes(
    const vector<IndexedMassSpectralPeak*>& xic,
    const Identification& identification, const int chargeState,
    unordered_map<string, vector<pair<double, double>>>&
        modifiedSequenceToIsotopicDistribution,
    const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks);

bool checkIsotopicEnvelopeCorrelation(
    map<int, vector<IsotopePeak>>& massShiftToIsotopePeaks,
    const IndexedMassSpectralPeak* peak, int chargeState,
    PpmTolerance& isotopeTolerance,
    const map<int, map<int, IndexedMassSpectralPeak>>& indexedPeaks);

struct filterResults {
    vector<double> expIntensity;
    vector<double> theorIntensity;
};

filterResults filterMassShiftToIsotopePeaks(vector<IsotopePeak>& isotopePeaks);

void setToNegativeOneIfNaN(double& value);

map<int, PSM> create_psm_map(const string& psm_file,
                             const string& psm_file_format = "percolator");

void cutPeak(ChromatographicPeak& peak, double identificationTime,
             unordered_map<string, vector<Ms1ScanInfo>>& _ms1Scans);

void runErrorChecking(const string& spectraFile, CruxLFQResults& lfqResults);

// void quantifyMatchBetweenRunsPeaks(const string& spectraFile, CruxLFQResults&
// lfqResults);

inline std::ostream& operator<<(std::ostream& os, const DetectionType& dt) {
    switch (dt) {
        case DetectionType::MSMS:
            os << "MSMS";
            break;
        case DetectionType::MBR:
            os << "MBR";
            break;
        case DetectionType::NotDetected:
            os << "NotDetected";
            break;
        case DetectionType::MSMSAmbiguousPeakfinding:
            os << "MSMSAmbiguousPeakfinding";
            break;
        case DetectionType::MSMSIdentifiedButNotQuantified:
            os << "MSMSIdentifiedButNotQuantified";
            break;
        case DetectionType::Imputed:
            os << "Imputed";
            break;
    }
    return os;
}

}  // namespace CruxLFQ

namespace std {
template <>
struct hash<CruxLFQ::Identification> {
    size_t operator()(const CruxLFQ::Identification& id) const {
        return hash<string>()(id.sequence) ^ hash<double>()(id.peptideMass) ^
               hash<double>()(id.monoIsotopicMass) ^
               hash<double>()(id.peakFindingMass) ^
               hash<int>()(id.precursorCharge) ^ hash<string>()(id.spectralFile) ^
               hash<double>()(id.ms2RetentionTimeInMinutes) ^
               hash<int>()(id.scanId) ^ hash<string>()(id.modifications);
    }
};

template <>
struct equal_to<CruxLFQ::Identification> {
    bool operator()(const CruxLFQ::Identification& id1,
                    const CruxLFQ::Identification& id2) const {
        // Compare id1 and id2 for equality
        // This uses the operator== that is already defined
        return id1 == id2;
    }
};
}  // namespace std