#include "CruxLFQApplication.h"

#include <cmath>
#include <exception>
#include <sstream>

#include "IndexedMassSpectralPeak.h"
#include "crux-lfq/IntensityNormalizationEngine.h"
#include "crux-lfq/Results.h"
#include "crux-lfq/Utils.h"
#include "io/carp.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/crux-utils.h"

using std::make_pair;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

using namespace pwiz::cv;
using namespace pwiz::msdata;
using pwiz::cv::MS_ms_level;
using pwiz::cv::MS_scan_start_time;
using pwiz::msdata::BinaryDataArrayPtr;
using pwiz::msdata::MSDataFile;
using pwiz::msdata::SpectrumListSimple;
using pwiz::msdata::SpectrumListSimplePtr;
using pwiz::msdata::SpectrumPtr;

typedef pwiz::msdata::SpectrumListPtr SpectrumListPtr;

int CruxLFQ::NUM_ISOTOPES_REQUIRED = 2;                       // Default value is 2
double CruxLFQ::PEAK_FINDING_PPM_TOLERANCE = 20.0;            // Default value is 20.0
double CruxLFQ::PPM_TOLERANCE = 10.0;                         // Default value is 10.0
bool CruxLFQ::ID_SPECIFIC_CHARGE_STATE = false;               // Default value is false
int CruxLFQ::MISSED_SCANS_ALLOWED = 1;                        // Default value is 1
double CruxLFQ::ISOTOPE_TOLERANCE_PPM = 5.0;                  // Default value is 5.0
bool CruxLFQ::INTEGRATE = false;                              // Default value is false
double CruxLFQ::DISCRIMINATION_FACTOR_TO_CUT_PEAK = 0.6;      // Default value is 0.6
bool CruxLFQ::QUANTIFY_AMBIGUOUS_PEPTIDES = false;            // Default value is false
bool CruxLFQ::USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT = false;  // Default value is false
bool CruxLFQ::NORMALIZE = false;                              // Default value is false

CruxLFQApplication::CruxLFQApplication() {}

CruxLFQApplication::~CruxLFQApplication() {}

int CruxLFQApplication::main(int argc, char** argv) {
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> spec_files = Params::GetStrings("spectrum files");
    return main(psm_file, spec_files);
}

int CruxLFQApplication::main(const string& psm_file, const vector<string>& spec_files) {
    carp(CARP_INFO, "Running crux-lfq...");

    CruxLFQ::NUM_ISOTOPES_REQUIRED = Params::GetInt("num-isotopes-required");                                   // Default value is 2
    CruxLFQ::PEAK_FINDING_PPM_TOLERANCE = Params::GetDouble("peak-finding-ppm-tolerance");                      // Default value is 20.0
    CruxLFQ::PPM_TOLERANCE = Params::GetDouble("ppm-tolerance");                                                // Default value is 10.0
    CruxLFQ::ID_SPECIFIC_CHARGE_STATE = Params::GetBool("id-specific-charge-state");                            // Default value is false
    CruxLFQ::MISSED_SCANS_ALLOWED = Params::GetInt("missed-scans-allowed");                                     // Default value is 1
    CruxLFQ::ISOTOPE_TOLERANCE_PPM = Params::GetDouble("isotope-tolerance-ppm");                                // Default value is 5.0
    CruxLFQ::INTEGRATE = Params::GetBool("integrate");                                                          // Default value is false
    CruxLFQ::DISCRIMINATION_FACTOR_TO_CUT_PEAK = Params::GetDouble("discrimination-factor-to-cut-peak");        // Default value is 0.6
    CruxLFQ::QUANTIFY_AMBIGUOUS_PEPTIDES = Params::GetBool("quantify-ambiguous-peptides");                      // Default value is false
    CruxLFQ::USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT = Params::GetBool("use-shared-peptides-for-protein-quant");  // Default value is false
    CruxLFQ::NORMALIZE = Params::GetBool("normalize");                                                          // Default value is false

    string output_dir = Params::GetString("output-dir");
    string psm_file_format = Params::GetString("psm-file-format");

    if (!FileUtils::Exists(psm_file)) {
        carp(CARP_FATAL, "PSM file %s not found", psm_file.c_str());
    }
    vector<CruxLFQ::PSM> psm_data = CruxLFQ::create_psm(psm_file, psm_file_format);
    CruxLFQ::CruxLFQResults lfqResults(spec_files);

    vector<CruxLFQ::Identification> allIdentifications;
    std::unordered_set<CruxLFQ::Identification> uniqueIdentifications;
    for (const string& spectra_file : spec_files) {
        vector<CruxLFQ::Identification> tempIdentifications = createIdentifications(psm_data, spectra_file);
        for (auto& id : tempIdentifications) {
            uniqueIdentifications.insert(id);
        }
    }
    std::copy(uniqueIdentifications.begin(), uniqueIdentifications.end(), std::back_inserter(allIdentifications));

    lfqResults.setPeptideModifiedSequencesAndProteinGroups(allIdentifications);

    unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxLFQ::calculateTheoreticalIsotopeDistributions(allIdentifications);

    vector<int> chargeStates = CruxLFQ::createChargeStates(allIdentifications);
    for (const string& spectra_file : spec_files) {
        SpectrumListPtr spectra_ms1 = loadSpectra(spectra_file, 1);
        carp(CARP_INFO, "Read %d spectra. for MS1. from %s", spectra_ms1->size(), spectra_file.c_str());

        CruxLFQ::IndexedSpectralResults indexResults = indexedMassSpectralPeaks(spectra_ms1, spectra_file);

        carp(CARP_INFO, "Finished indexing peaks for %s", spectra_file.c_str());

        // TODO Continue from this function
        CruxLFQ::quantifyMs2IdentifiedPeptides(
            spectra_file,
            allIdentifications,
            chargeStates,
            indexResults._ms1Scans,
            indexResults._indexedPeaks,
            modifiedSequenceToIsotopicDistribution,
            lfqResults);
        CruxLFQ::runErrorChecking(spectra_file, lfqResults);

        carp(CARP_INFO, "Finished processing %s", spectra_file.c_str());
    }

    if (CruxLFQ::NORMALIZE) {
        CruxLFQ::IntensityNormalizationEngine intensityNormalizationEngine(
            lfqResults,
            CruxLFQ::INTEGRATE,
            CruxLFQ::QUANTIFY_AMBIGUOUS_PEPTIDES);
        intensityNormalizationEngine.NormalizeResults();
    }

    lfqResults.calculatePeptideResults(CruxLFQ::QUANTIFY_AMBIGUOUS_PEPTIDES);
    lfqResults.calculateProteinResultsMedianPolish(CruxLFQ::USE_SHARED_PEPTIDES_FOR_PROTEIN_QUANT);
    const std::string mod_pep_results_file = make_file_path("crux-lfq-mod-pep.txt");
    const std::string peak_results_file = make_file_path("crux-lfq-peaks.txt");
    lfqResults.writeResults(mod_pep_results_file, peak_results_file, spec_files);

    return 0;
}

string CruxLFQApplication::getName() const {
    return "crux-lfq";
}

string CruxLFQApplication::getDescription() const {
    return "[[nohtml:This command reads a set of PSMs and a corresponding set of spectrum files"
           "and carries out label-free quantification (LFQ) for each detected peptide.]]"
           "[[html:<p>This command reads a set of PSMs and a corresponding set of spectrum files "
           "and carries out label-free quantification (LFQ) for each detected peptide."
           "The algorithm follows that of FlashLFQ: "
           "Millikin RJ, Solntsev SK, Shortreed MR, Smith LM. &quot;<a href=\""
           "https://pubmed.ncbi.nlm.nih.gov/29083185/\">Ultrafast Peptide Label-Free Quantification with FlashLFQ.</a>&quot;"
           "<em>Journal of Proteome Research</em>. 17(1):386-391, 2018.</blockquote><p>]]";
}

vector<string> CruxLFQApplication::getArgs() const {
    string arr[] = {
        "lfq-peptide-spectrum matches",
        "spectrum files+"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> CruxLFQApplication::getOptions() const {
    string arr[] = {
        "score",
        "threshold",
        "smaller-is-better",
        "fileroot",
        "output-dir",
        "overwrite",
        "parameter-file",
        "verbosity",
        "num-isotopes-required",
        "peak-finding-ppm-tolerance",
        "ppm-tolerance",
        "id-specific-charge-state",
        "missed-scans-allowed",
        "isotope-tolerance-ppm",
        "integrate",
        "discrimination-factor-to-cut-peak",
        "quantify-ambiguous-peptides",
        "use-shared-peptides-for-protein-quant",
        "normalize",
    };
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<pair<string, string>> CruxLFQApplication::getOutputs() const {
    vector<pair<string, string>> outputs;
    outputs.push_back(make_pair("crux-lfq.txt",
                                "A tab-delimited text file in which rows are peptides, "
                                "columns correspond to the different spectrum files, "
                                "and values are peptide quantifications.  "
                                "If a peptide is not detected in a given run, "
                                "then its corresponding quantification value is NaN."));
    outputs.push_back(make_pair("crux-lfq.params.txt",
                                "A file containing the name and value of all parameters/options"
                                " for the current operation. Not all parameters in the file may have"
                                " been used in the operation. The resulting file can be used with the "
                                "--parameter-file option for other Crux programs."));
    outputs.push_back(make_pair("crux-lfq.log.txt",
                                "A log file containing a copy of all messages that were printed to the screen during execution."));
    return outputs;
}

COMMAND_T CruxLFQApplication::getCommand() const {
    return CRUX_LFQ_COMMAND;
}

bool CruxLFQApplication::needsOutputDirectory() const {
    return true;
}

// TODO: Add parameter processing
void CruxLFQApplication::processParams() {
}

SpectrumListPtr CruxLFQApplication::loadSpectra(const string& file, int msLevel) {
    try {
        MSDataFile msd(file);
        SpectrumListPtr originalSpectrumList = msd.run.spectrumListPtr;
        if (!originalSpectrumList) {
            carp(CARP_FATAL, "Error reading spectrum file %s", file.c_str());
        }
        SpectrumListSimplePtr filteredSpectrumList(new SpectrumListSimple);

        for (size_t i = 0; i < originalSpectrumList->size(); ++i) {
            const bool getBinaryData = true;
            SpectrumPtr spectrum = originalSpectrumList->spectrum(i, getBinaryData);

            int spectrumMSLevel = spectrum->cvParam(MS_ms_level).valueAs<int>();
            if (spectrumMSLevel == msLevel) {
                // Add the spectrum to the filtered list
                filteredSpectrumList->spectra.push_back(spectrum);
            }
        }

        return filteredSpectrumList;
    } catch (const std::exception& e) {
        carp(CARP_INFO, "Error:  %s", e.what());
        return nullptr;
    }
}

IndexedSpectralResults CruxLFQApplication::indexedMassSpectralPeaks(SpectrumListPtr spectrum_collection, const string& spectra_file) {
    string _spectra_file(spectra_file);

    vector<vector<IndexedMassSpectralPeak>> _indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;

    _ms1Scans[_spectra_file] = vector<Ms1ScanInfo>();

    IndexedSpectralResults index_results{_indexedPeaks, _ms1Scans};

    if (!spectrum_collection) {
        return index_results;
    }

    double maxMz = 0.0;
    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        const bool getBinaryData = true;
        SpectrumPtr spectrum = spectrum_collection->spectrum(i, getBinaryData);
        if (spectrum) {
            BinaryDataArrayPtr mzs = spectrum->getMZArray();
            if (mzs) {
                const std::vector<double>& mzArray = mzs->data;
                for (size_t j = 0; j < mzArray.size(); ++j) {
                    double mz = mzArray[j];
                    if (mz > maxMz) {
                        maxMz = mz;
                    }
                }
            }
        }
    }
    int size = (int)std::ceil(maxMz * BINS_PER_DALTON) + 1;
    index_results._indexedPeaks.resize(size);

    int scanIndex = 0;
    int oneBasedScanNumber = 1;

    for (size_t i = 0; i < spectrum_collection->size(); ++i) {
        const bool getBinaryData = true;
        SpectrumPtr spectrum = spectrum_collection->spectrum(i, getBinaryData);

        if (spectrum) {
            BinaryDataArrayPtr mzs = spectrum->getMZArray();
            BinaryDataArrayPtr intensities = spectrum->getIntensityArray();

            double retentionTime = spectrum->scanList.scans[0].cvParam(MS_scan_start_time).timeInSeconds();

            if (mzs && intensities) {
                const std::vector<double>& mzArray = mzs->data;
                const std::vector<double>& intensityArray = intensities->data;
                for (size_t j = 0; j < mzArray.size(); ++j) {
                    double mz = mzArray[j];
                    int roundedMz = static_cast<int>(std::round(mz * BINS_PER_DALTON));
                    IndexedMassSpectralPeak spec_data(
                        mz,                 // mz value
                        intensityArray[j],  // intensity
                        scanIndex,          // zeroBasedMs1ScanIndex
                        retentionTime);

                    if (index_results._indexedPeaks[roundedMz].empty()) {
                        index_results._indexedPeaks[roundedMz] = vector<IndexedMassSpectralPeak>();
                    }
                    index_results._indexedPeaks[roundedMz].push_back(spec_data);
                    Ms1ScanInfo scan = {oneBasedScanNumber, scanIndex, retentionTime};
                    index_results._ms1Scans[spectra_file].push_back(scan);
                }
            }

            scanIndex++;
            oneBasedScanNumber++;
        }
    }
    return index_results;
}

// Make this a multithreaded process
vector<Identification> CruxLFQApplication::createIdentifications(const vector<PSM>& psm_data, const string& spectra_file) {
    carp(CARP_INFO, "Creating indentifications, this may take a bit of time, do not terminate the process...");

    vector<Identification> allIdentifications;
    string _spectra_file(spectra_file);

    for (const auto psm : psm_data) {
        Identification identification;

        identification.sequence = psm.sequence_col;
        identification.monoIsotopicMass = psm.monoisotopic_mass_col;
        identification.peptideMass = psm.peptide_mass_col;
        identification.precursorCharge = psm.charge_col;
        identification.spectralFile = _spectra_file;
        identification.ms2RetentionTimeInMinutes = psm.retention_time;
        identification.scanId = psm.scan_col;
        identification.modifications = psm.modifications;
        allIdentifications.push_back(identification);
    }

    return allIdentifications;
}