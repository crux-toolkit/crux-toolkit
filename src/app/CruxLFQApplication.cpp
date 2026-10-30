#include "CruxLFQApplication.h"

#include <cmath>
#include <exception>
#include <fstream>
#include <list>
#include <memory>
#include <sstream>

#include "IndexedMassSpectralPeak.h"
#include "LFQMetaData.h"
#include "app/tide/mass_constants.h"
#include "app/tide/modifications.h"
#include "crux-lfq/IntensityNormalizationEngine.h"
#include "crux-lfq/Results.h"
#include "crux-lfq/Utils.h"
#include "io/carp.h"
#include "model/Peptide.h"
#include "util/FileUtils.h"
#include "util/Params.h"
#include "util/crux-utils.h"

using std::list;
using std::make_pair;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

namespace {
    string fileNameWithoutExtension(const string& path) {
        size_t lastSlash = path.find_last_of("/\\");
        string name = (lastSlash == string::npos) ? path : path.substr(lastSlash + 1);
        size_t lastDot = name.find_last_of('.');
        return (lastDot == string::npos) ? name : name.substr(0, lastDot);
    }
}

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
int CruxLFQ::MaxThreads = 1;                                  // Default value is 1

CruxLFQApplication::CruxLFQApplication() {}

CruxLFQApplication::~CruxLFQApplication() {}

int CruxLFQApplication::main(int argc, char** argv) {
    string psm_file = Params::GetString("lfq-peptide-spectrum matches");
    vector<string> spec_files = Params::GetStrings("spectrum files");
    string specfile_replicates = Params::GetString("specfile-replicates");
    return main(psm_file, spec_files, specfile_replicates);
}

int CruxLFQApplication::main(const string& psm_file, const vector<string>& spec_files, const string& specfile_replicates) {
    return main(psm_file, spec_files, specfile_replicates, unordered_map<string, string>());
}

int CruxLFQApplication::main(const string& psm_file, const vector<string>& spec_files, const string& specfile_replicates,
                             const unordered_map<string, string>& file_aliases) {
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
    CruxLFQ::MaxThreads = Params::GetInt("num-threads");                                                        // Default value is 1

    string output_dir = Params::GetString("output-dir");

    if (!FileUtils::Exists(psm_file)) {
        carp(CARP_FATAL, "PSM file %s not found", psm_file.c_str());
    }

    string psm_file_format = Params::GetString("psm-file-format");
    vector<CruxLFQ::PSM> psm_data;
    if (psm_file_format == "percolator") {
        psm_data = create_percolator_psm(psm_file);
    } else {
        psm_data = CruxLFQ::create_psm(psm_file);
    }

    if (CruxLFQ::NORMALIZE && !FileUtils::Exists(specfile_replicates)) {
        carp(CARP_FATAL, "Normalization requires a spectrum file replicates file. File '%s' not found.", specfile_replicates.c_str());
    }
    CruxLFQ::CruxLFQResults lfqResults = CruxLFQ::NORMALIZE
        ? CruxLFQ::CruxLFQResults(specfile_replicates)
        : CruxLFQ::CruxLFQResults(spec_files);

    vector<CruxLFQ::Identification> allIdentifications;
    std::unordered_set<CruxLFQ::Identification> uniqueIdentifications;
    vector<CruxLFQ::Identification> tempIdentifications = createIdentifications(psm_data, spec_files, file_aliases);
    for (auto& id : tempIdentifications) {
        uniqueIdentifications.insert(id);
    }
    std::copy(uniqueIdentifications.begin(), uniqueIdentifications.end(), std::back_inserter(allIdentifications));

    lfqResults.setPeptideModifiedSequencesAndProteinGroups(allIdentifications);

    unordered_map<string, vector<pair<double, double>>> modifiedSequenceToIsotopicDistribution = CruxLFQ::calculateTheoreticalIsotopeDistributions(allIdentifications);

    vector<int> chargeStates = CruxLFQ::createChargeStates(allIdentifications);

    // {
    //     const std::string id_file = make_file_path("crux-lfq-identifications.txt");
    //     std::ofstream id_out(id_file);
    //     if (id_out.is_open()) {
    //         id_out << "sequence\tpeptideMass\tmonoIsotopicMass\tpeakFindingMass\t"
    //                   "precursorCharge\tspectralFile\tms2RetentionTimeInMinutes\t"
    //                   "scanId\tmodifications\tprotein_id\n";
    //         for (const auto& id : allIdentifications) {
    //             id_out << id.sequence << "\t"
    //                    << id.peptideMass << "\t"
    //                    << id.monoIsotopicMass << "\t"
    //                    << id.peakFindingMass << "\t"
    //                    << id.precursorCharge << "\t"
    //                    << id.spectralFile << "\t"
    //                    << id.ms2RetentionTimeInMinutes << "\t"
    //                    << id.scanId << "\t"
    //                    << id.modifications << "\t"
    //                    << id.protein_id << "\n";
    //         }
    //     } else {
    //         carp(CARP_WARNING, "Could not open identifications file for writing: %s", id_file.c_str());
    //     }
    //     // Exit after writing identifications file if no spectrum files are provided, since there's nothing more to do
    //     return 0;
    // }

    for (const string& spectra_file : spec_files) {
        Crux::SpectrumCollection* spectra_ms1 = loadSpectra(spectra_file, 1);
        carp(CARP_INFO, "Read %d spectra. for MS1. from %s", spectra_ms1->getNumSpectra(), spectra_file.c_str());
        if (spectra_ms1->getNumSpectra() == 0) {
            carp(CARP_FATAL, "No MS1 spectra found in '%s'. LFQ quantifies from MS1 scans, so it needs a "
                 "spectrum file that contains them (e.g. mzML or mzXML, not ms2 or mgf).", spectra_file.c_str());
        }

        CruxLFQ::IndexedSpectralResults indexResults = indexedMassSpectralPeaks(spectra_ms1, spectra_file);

        carp(CARP_INFO, "Finished indexing peaks for %s", spectra_file.c_str());

        vector<CruxLFQ::Identification> filteredIdentifications;

        std::copy_if(
            allIdentifications.begin(),
            allIdentifications.end(),
            std::back_inserter(filteredIdentifications),
            [&spectra_file](const Identification& identification) {
                return identification.spectralFile == spectra_file;
            });
        carp(CARP_INFO, "Filtered %zu identifications out of %zu total for '%s'",
             filteredIdentifications.size(), allIdentifications.size(), spectra_file.c_str());

        vector<vector<CruxLFQ::IndexedMassSpectralPeak>*> convertedPeaks;

        auto metadata = &CruxLFQ::LFQMetaData::getInstance();
        std::transform(
            indexResults._indexedPeaks.begin(),
            indexResults._indexedPeaks.end(),
            std::back_inserter(convertedPeaks),
            [](vector<CruxLFQ::IndexedMassSpectralPeak>& innerVec) {
                return &innerVec;
            });
        metadata->setIndexedPeaks(&convertedPeaks);
        metadata->setMs1Scans(&indexResults._ms1Scans);

        CruxLFQ::quantifyMs2IdentifiedPeptides(
            spectra_file,
            filteredIdentifications,
            chargeStates,
            modifiedSequenceToIsotopicDistribution,
            &lfqResults);
        CruxLFQ::runErrorChecking(spectra_file, lfqResults);

        carp(CARP_INFO, "Finished processing %s", spectra_file.c_str());
        delete spectra_ms1;
    }

    if (CruxLFQ::NORMALIZE) {
        CruxLFQ::IntensityNormalizationEngine intensityNormalizationEngine(
            lfqResults,
            CruxLFQ::INTEGRATE,
            false,
            output_dir);
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
    return "lfq";
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
        "specfile-replicates",
        "psm-file-format",
        "is-rt-seconds",
        "spectrum-parser",
        "num-threads",
        "mods-spec",
        "nterm-peptide-mods-spec",
        // "nterm-protein-mods-spec",
        "cterm-peptide-mods-spec",
        // "cterm-protein-mods-spec",
        "lfq-q-value-threshold",
        "is-psm-filtered"};
    return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<pair<string, string>> CruxLFQApplication::getOutputs() const {
    vector<pair<string, string>> outputs;
    outputs.push_back(make_pair("crux-lfq-mod-pep.txt",
                                "A tab-delimited text file in which rows are peptides, "
                                "columns correspond to the different spectrum files, "
                                "and values are peptide quantifications.  "
                                "If a peptide is not detected in a given run, "
                                "then its corresponding quantification value is NaN."));
    outputs.push_back(make_pair("crux-lfq-peaks.txt",
                                "A tab-delimited text file in which rows are peaks, "
                                "columns correspond to meta-data about the peaks"));
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

Crux::SpectrumCollection* CruxLFQApplication::loadSpectra(const string& file, int msLevel) {
    Crux::SpectrumCollection* spectra(SpectrumCollectionFactory::create(file.c_str()));
    spectra->parse(msLevel);
    return spectra;
}

IndexedSpectralResults CruxLFQApplication::indexedMassSpectralPeaks(Crux::SpectrumCollection* spectrum_collection, const string& spectra_file) {
    // string _spectra_file(spectra_file);

    vector<vector<IndexedMassSpectralPeak>> _indexedPeaks;
    unordered_map<string, vector<Ms1ScanInfo>> _ms1Scans;

    _ms1Scans[spectra_file] = vector<Ms1ScanInfo>();

    IndexedSpectralResults index_results{_indexedPeaks, _ms1Scans};

    if (spectrum_collection->getNumSpectra() <= 0) {
        return index_results;
    }

    int scanIndex = 0;
    for (auto spectrum = spectrum_collection->begin(); spectrum != spectrum_collection->end(); ++spectrum) {
        if (*spectrum != nullptr) {
            double retentionTime = (*spectrum)->getRTime();
            int nativeScanNumber = (*spectrum)->getFirstScan();  // Native scan number from mzML
            for (auto peak = (*spectrum)->begin(); peak != (*spectrum)->end(); ++peak) {
                if (*peak != nullptr) {
                    double mz = (*peak)->getLocation();
                    int roundedMz = static_cast<int>(std::round(mz * BINS_PER_DALTON));

                    if (index_results._indexedPeaks.size() <= roundedMz) {
                        int size = roundedMz + 1;
                        index_results._indexedPeaks.resize(size);
                    }

                    index_results._indexedPeaks[roundedMz].emplace_back(
                        mz,                       // mz value
                        (*peak)->getIntensity(),  // intensity
                        scanIndex,                // zeroBasedMs1ScanIndex
                        nativeScanNumber,         // nativeScanNumber from mzML
                        retentionTime             // retentionTime
                    );
                }
            }
            index_results._ms1Scans[spectra_file].emplace_back(nativeScanNumber, scanIndex, retentionTime);
            scanIndex++;
        }
    }

    return index_results;
}

// Make this a multithreaded process
vector<Identification> CruxLFQApplication::createIdentifications(const vector<PSM>& psm_data, const vector<string>& spec_files,
                                                                 const unordered_map<string, string>& file_aliases) {
    carp(CARP_INFO, "Creating identifications, this may take a bit of time, do not terminate the process...");

    // Build a map from file stem -> CLI spectra file path (like C# FindMatchingSpecFile)
    unordered_map<string, string> stemToSpecFile;
    for (const auto& spec_file : spec_files) {
        stemToSpecFile[fileNameWithoutExtension(spec_file)] = spec_file;
    }

    // Resolve a PSM's file entry to a spectra file: directly by stem, or via an
    // alias (e.g. the search results file a pipeline produced from that spectra file).
    auto resolveSpecFile = [&](const string& file_name) -> const string* {
        auto it = stemToSpecFile.find(fileNameWithoutExtension(file_name));
        if (it != stemToSpecFile.end()) {
            return &it->second;
        }
        auto alias = file_aliases.find(file_name);
        if (alias == file_aliases.end()) {
            alias = file_aliases.find(FileUtils::BaseName(file_name));
        }
        if (alias != file_aliases.end()) {
            it = stemToSpecFile.find(fileNameWithoutExtension(alias->second));
            if (it != stemToSpecFile.end()) {
                return &it->second;
            }
        }
        return NULL;
    };

    // MS2 retention times by scan, loaded only for files whose PSMs lack RT.
    unordered_map<string, unordered_map<int, double>> ms2RetentionTimes;
    auto ms2RetentionTime = [&](const string& spec_file, int scan, double* rt) -> bool {
        auto file_it = ms2RetentionTimes.find(spec_file);
        if (file_it == ms2RetentionTimes.end()) {
            carp(CARP_INFO, "Reading MS2 retention times from '%s' for PSMs that lack them", spec_file.c_str());
            Crux::SpectrumCollection* spectra_ms2 = loadSpectra(spec_file, 2);
            unordered_map<int, double>& scanToRt = ms2RetentionTimes[spec_file];
            for (auto spectrum = spectra_ms2->begin(); spectrum != spectra_ms2->end(); ++spectrum) {
                scanToRt[(*spectrum)->getFirstScan()] = (*spectrum)->getRTime();
            }
            delete spectra_ms2;
            file_it = ms2RetentionTimes.find(spec_file);
        }
        auto scan_it = file_it->second.find(scan);
        if (scan_it == file_it->second.end()) {
            return false;
        }
        *rt = scan_it->second;
        return true;
    };

    vector<Identification> allIdentifications;
    size_t unmatched = 0;

    for (const auto& psm : psm_data) {
        const string* specFile = resolveSpecFile(psm.file_name);
        if (specFile == NULL) {
            // PSM is from a file we're not processing — skip it
            ++unmatched;
            continue;
        }

        // A missing retention time (absent, NaN or non-positive in the PSM file)
        // is taken from the MS2 scan itself, so it is in the same units as the
        // MS1 scans it will be matched against.
        double retentionTime = psm.retention_time;
        if (!(retentionTime > 0) && !ms2RetentionTime(*specFile, psm.scan_col, &retentionTime)) {
            carp(CARP_FATAL, "PSM for scan %d has no retention time, and scan %d was not found in '%s'.",
                 psm.scan_col, psm.scan_col, specFile->c_str());
        }

        allIdentifications.emplace_back(
            psm.sequence_col,
            psm.monoisotopic_mass_col,
            psm.peptide_mass_col,
            psm.charge_col,
            *specFile,  // resolved CLI spectra file path
            retentionTime,
            psm.scan_col,
            psm.modifications,
            psm.protein_id);
    }

    carp(CARP_INFO, "Created %zu identifications matching %zu spectra files (from %zu total PSMs)",
         allIdentifications.size(), spec_files.size(), psm_data.size());
    if (!psm_data.empty() && allIdentifications.empty()) {
        carp(CARP_FATAL, "None of the %zu PSMs could be matched to the spectrum files given. The PSM file's "
             "file column must name the spectrum files (e.g. '%s').",
             psm_data.size(), spec_files.empty() ? "" : spec_files.front().c_str());
    } else if (unmatched > 0) {
        carp(CARP_WARNING, "%zu of %zu PSMs are from spectrum files that were not given and will not be quantified.",
             unmatched, psm_data.size());
    }

    return allIdentifications;
}

void CruxLFQApplication::gen_mods(
    string sequence_col,
    ModPosition mod_psn_type,
    const pb::ModTable* mod_table_,
    vector<Crux::Modification>& mods) {
    double mod_mass;
    ModPosition position = mod_psn_type;

    for (size_t i = 0; i < sequence_col.size(); ++i) {
        char AA = sequence_col[i];
        for (int i = 0; i < mod_table_->static_mod_size(); i++) {
            const pb::Modification& mod = mod_table_->static_mod(i);
            const ModificationDefinition* mod_;
            if (mod.has_delta() && mod.has_amino_acids() && mod.has_name()) {
                string AAs = mod.amino_acids();
                int AA_len = AAs.length();
                for (int j = 0; j < AA_len; ++j) {
                    if (AAs[j] == AA || AAs[j] == 'X') {  // Found a static mod for Amino acid AA;
                        mod_mass = mod.delta();
                        mod_ = ModificationDefinition::New(string(1, AAs[j]), mod_mass, position, true);
                        mods.push_back(Crux::Modification(mod_, i));
                    }
                }
            }
        }
    }
}

vector<PSM> CruxLFQApplication::create_percolator_psm(const string& psm_file) {
    bool is_rt_seconds = Params::GetBool("is-rt-seconds");
    double q_value_threshold = Params::GetDouble("lfq-q-value-threshold");
    bool filtered = Params::GetBool("is-psm-filtered");

    string mods_spec = Params::GetString("mods-spec");
    // if (!mods_spec.empty()) {
    //     if (std::isdigit(mods_spec[0])) {
    //         carp(CARP_FATAL, "mods-spec must be static not variable");
    //     }
    // } else {
    //     carp(CARP_FATAL, "mods-spec can't be empty for percolator PSM file formats");
    // }

    if (mods_spec.empty()) {
        carp(CARP_FATAL, "mods-spec can't be empty for percolator PSM file formats");
    }

    vector<PSM> psm_data;
    std::ifstream file(psm_file);
    if (!file.is_open()) {
        carp(CARP_FATAL, "Error: Could not open the PSM file!");
    }

    string sequence_col, protein_id, line, psm_id, file_name_col;
    int scan_col, charge_col;
    double peptide_mass_col, q_value, retention_time;

    // Read the header line and ignore
    std::getline(file, line);

    VariableModTable var_mod_table;
    var_mod_table.ClearTables();

    // parse regular amino acid modifications
    if (!var_mod_table.Parse(mods_spec.c_str())) {
        carp(CARP_FATAL, "Error parsing mods");
    }
    // parse terminal modifications
    mods_spec = Params::GetString("cterm-peptide-mods-spec");
    if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPEP)) {
        carp(CARP_FATAL, "Error parsing c-terminal peptide mods");
    }

    mods_spec = Params::GetString("nterm-peptide-mods-spec");
    if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPEP)) {
        carp(CARP_FATAL, "Error parsing n-terminal peptide mods");
    }

    // mods_spec = Params::GetString("cterm-protein-mods-spec");
    // if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), CTPRO)) {
    //     carp(CARP_FATAL, "Error parsing c-terminal protein mods");
    // }

    // mods_spec = Params::GetString("nterm-protein-mods-spec");
    // if (!mods_spec.empty() && !var_mod_table.Parse(mods_spec.c_str(), NTPRO)) {
    //     carp(CARP_FATAL, "Error parsing n-terminal protein mods");
    // }
    var_mod_table.SerializeUniqueDeltas();
    if (!MassConstants::Init(var_mod_table.ParsedModTable(),
                             var_mod_table.ParsedNtpepModTable(),
                             var_mod_table.ParsedCtpepModTable(),
                             var_mod_table.ParsedNtproModTable(),
                             var_mod_table.ParsedCtproModTable(),
                             MassConstants::bin_width_,
                             MassConstants::bin_offset_)) {
        carp(CARP_FATAL, "Error in MassConstants::Init");
    }

    // The header was consumed above; line numbers below count it as line 1.
    int line_number = 1;
    while (std::getline(file, line)) {
        line_number++;
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;

        // Split line by tabs
        while (std::getline(iss, token, '\t')) {
            tokens.push_back(token);
        }

        // Debug output for the failing line
        if (tokens.size() != 8) {
            carp(CARP_ERROR, "Line %d has %zu tokens (expected 8): [%s]",
                 line_number, tokens.size(), line.c_str());

            // Print each token with its length
            for (size_t i = 0; i < tokens.size(); ++i) {
                carp(CARP_ERROR, "  Token %zu (len=%zu): [%s]",
                     i, tokens[i].length(), tokens[i].c_str());
            }
        }

        if (tokens.size() < 8) {
            carp(CARP_FATAL, "PSM file has malformed data on line %d: %s",
                 line_number, line.c_str());
        }
        psm_id = tokens[0];
        file_name_col = tokens[1];
        q_value = std::stod(tokens[3]);
        sequence_col = tokens[5];
        protein_id = tokens[6];
        retention_time = std::stod(tokens.back());

        // const char* booleanText = filtered ? "true" : "false";

        // carp(CARP_INFO, "q_value %f, q_value_threshold %f, filtered %s", q_value, q_value_threshold, booleanText);
        if (!filtered && q_value > q_value_threshold) {
            // carp(CARP_INFO, "I was passed");
            continue;
        }
        if (is_rt_seconds) {
            retention_time = retention_time / 60.0;
        }

        std::istringstream ss(psm_id);
        std::string psm_id_tokens;
        std::vector<std::string> parts;

        while (std::getline(ss, psm_id_tokens, '_')) {
            parts.push_back(psm_id_tokens);
        }

        if (parts.size() == 5) {
            scan_col = std::stoi(parts[2]);
            charge_col = std::stoi(parts[3]);
        } else {
            carp(CARP_FATAL, "Error: Unexpected PSMId format.");
        }

        vector<Crux::Modification> mods;
        Crux::Modification::FromSeq(sequence_col, NULL, &mods);

        if (sequence_col.size() > 4) {  // Ensure the string is long enough
            sequence_col = sequence_col.substr(2, sequence_col.size() - 4);
        }

        gen_mods(sequence_col, ANY, MassConstants::mod_table_, mods);
        gen_mods(sequence_col, PEPTIDE_N, MassConstants::n_mod_table_, mods);
        gen_mods(sequence_col, PEPTIDE_C, MassConstants::c_mod_table_, mods);

        Crux::Peptide* peptide = new Crux::Peptide();
        string unmodSeq = Crux::Peptide::unmodifySequence(sequence_col);
        peptide->setUnmodifiedSequence(unmodSeq);
        peptide->setMods(mods);
        peptide_mass_col = peptide->calcModifiedMass();
        // peptide_mass_col = peptide->calcMass(MONO);
        // peptide_mass_col = peptide->calcMass(AVERAGE);

        // carp(CARP_INFO, "peptide_mass_col %f", peptide_mass_col);

        psm_data.emplace_back(sequence_col,
                              scan_col,
                              charge_col,
                              peptide_mass_col,
                              peptide_mass_col,
                              sequence_col,
                              retention_time,
                              protein_id,
                              file_name_col);
    }
    return psm_data;
}