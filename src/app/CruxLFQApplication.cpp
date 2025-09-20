#include "CruxLFQApplication.h"

#include <cmath>
#include <exception>
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
    string specfile_replicates = Params::GetString("spectrum file replicates");
    return main(psm_file, spec_files, specfile_replicates);
}

int CruxLFQApplication::main(const string& psm_file, const vector<string>& spec_files, const string& specfile_replicates) {
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

    CruxLFQ::CruxLFQResults lfqResults(spec_files);
    if (CruxLFQ::NORMALIZE && !FileUtils::Exists(specfile_replicates)) {
        carp(CARP_INFO, "Spectrum file replicates file %s not found", specfile_replicates.c_str());
        carp(CARP_FATAL, "Normalization requires a spectrum file replicates file.");
        lfqResults = CruxLFQ::CruxLFQResults(specfile_replicates);
    }

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
        Crux::SpectrumCollection* spectra_ms1 = loadSpectra(spectra_file, 1);
        carp(CARP_INFO, "Read %d spectra. for MS1. from %s", spectra_ms1->getNumSpectra(), spectra_file.c_str());

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
        "spectrum files+",
        "spectrum file replicates"};
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
    int oneBasedScanNumber = 1;
    for (auto spectrum = spectrum_collection->begin(); spectrum != spectrum_collection->end(); ++spectrum) {
        if (*spectrum != nullptr) {
            double retentionTime = (*spectrum)->getRTime();
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
                        retentionTime             // retentionTime
                    );
                }
            }
            index_results._ms1Scans[spectra_file].emplace_back(oneBasedScanNumber, scanIndex, retentionTime);
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

    for (const auto& psm : psm_data) {
        allIdentifications.emplace_back(
            psm.sequence_col,
            psm.monoisotopic_mass_col,
            psm.peptide_mass_col,
            psm.charge_col,
            spectra_file,
            psm.retention_time,
            psm.scan_col,
            psm.modifications,
            psm.protein_id);
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

    string sequence_col, protein_id, line, psm_id;
    int scan_col, charge_col;
    double peptide_mass_col, q_value, retention_time;

    // Read the header line and ignore
    std::getline(file, line);

    VariableModTable var_mod_table;
    var_mod_table.ClearTables();

    // parse regular amino acid modifications
    carp(CARP_DEBUG, "mods_spec='%s'", mods_spec.c_str());
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
    carp(CARP_INFO, "%s", mods_spec.c_str());
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

    int line_number = 0;
    while (std::getline(file, line)) {
        line_number++;
        // Skip header
        if (line_number == 1) continue;
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
                              protein_id);
    }
    return psm_data;
}