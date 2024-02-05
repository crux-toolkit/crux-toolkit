#include <iomanip>

#include "TideLiteMatchSet.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "tide/peptide_lite.h"
#include "crux_version.h"

// SCORE_FUNCTION_T is defined in ./src/model/objects.h
SCORE_FUNCTION_T TideLiteMatchSet::curScoreFunction_ = INVALID_SCORE_FUNCTION;
int TideLiteMatchSet::top_matches_ = 0;
int TideLiteMatchSet::decoy_num_ = 0;
int TideLiteMatchSet::mass_precision_ = 0;
int TideLiteMatchSet::score_precision_ = 0;
int TideLiteMatchSet::mod_precision_ = 0;
bool TideLiteMatchSet::concat_ = false;
string TideLiteMatchSet::decoy_prefix_ = "";
int TideLiteMatchSet::psm_id_mzTab_  = 1;
string TideLiteMatchSet::fasta_file_name_ = "null";


// column IDs are defined in ./src/io/MatchColumns.h and .cpp
// The order of the columns in the output is defined here by the order of the column Ids
int TideLiteMatchSet::XCorr_tsv_cols[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, RETENTION_TIME_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
    XCORR_RANK_COL, DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
 int TideLiteMatchSet::Pvalues_tsv_cols[] = {  //TODO: update the colums.
    FILE_COL, SCAN_COL, CHARGE_COL, RETENTION_TIME_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
    RESIDUE_EVIDENCE_COL, RESIDUE_PVALUE_COL, BOTH_PVALUE_COL, BOTH_PVALUE_RANK, 
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
int TideLiteMatchSet::Diameter_tsv_cols[] = { //TODO: update the colums.
    FILE_COL, SCAN_COL, CHARGE_COL, RETENTION_TIME_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, // TODO: add a variable for the matching theoretical peak series.
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    

int TideLiteMatchSet::XCorr_mzTab_cols[] = {
    MZTAB_PSH, MZTAB_SEQUENCE, MZTAB_PSM_ID, MZTAB_ACCESSION, MZTAB_UNIQUE, MZTAB_DATABASE,
    MZTAB_DATABASE_VERSION, MZTAB_SEARCH_ENGINE, 
    MZTAB_SEARCH_ENGINE_SCORE_1,  // [MS, MS:1001155, The SEQUEST result 'XCorr'.]
    MZTAB_SEARCH_ENGINE_SCORE_2,  // [MS, MS:1001143, The SEQUEST result 'DeltaCn'.]
    MZTAB_SEARCH_ENGINE_SCORE_3,  // [MS, MS:1003358, XCorr rank]
    MZTAB_SEARCH_ENGINE_SCORE_4,  // [MS, MS:1003359, exact p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_5,  // [MS, MS:1003360, refactored XCorr'.]
    MZTAB_SEARCH_ENGINE_SCORE_6,  // [MS, MS:1003361, res-ev score'.]
    MZTAB_SEARCH_ENGINE_SCORE_7,  // [MS, MS:1003362, res-ev rank'.]
    MZTAB_SEARCH_ENGINE_SCORE_8,  // [MS, MS:1003363, res-ev p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_9,  // [MS, MS:1003364, combined p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_10, // [MS, MS:1003365, combined p-value rank'.]
    MZTAB_SEARCH_ENGINE_SCORE_11, // [MS, MS:1003366, tailor score'.]
    MZTAB_MODIFICATIONS, MZTAB_RETENTION_TIME,
    MZTAB_CHARGE, MZTAB_EXP_MASS_TO_CHARGE, MZTAB_CALC_MASS_TO_CHARGE, MZTAB_SPECTRA_REF,
    MZTAB_PRE, MZTAB_POST, MZTAB_START, MZTAB_END, MZTAB_OPT_MS_RUN_1_SPECTRUM_NEUTRAL_MASS,
    MZTAB_OPT_MS_RUN_1_DELTA_LCN, MZTAB_OPT_MS_RUN_1_DISTINCT_MATCHES_PER_SPEC,
    MZTAB_OPT_MS_RUN_1_TARGET_DECOY, MZTAB_OPT_MS_RUN_1_ORIGINAL_TARGET_SEQUENCE_COL, 
    MZTAB_OPT_MS_RUN_1_DECOY_INDEX
  };    

int TideLiteMatchSet::XCorr_pin_cols[] = {
    POUT_PSMID_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
    XCORR_RANK_COL, DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
 int TideLiteMatchSet::Pvalues_pin_cols[] = {  //TODO: update the colums.
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
    RESIDUE_EVIDENCE_COL, RESIDUE_PVALUE_COL, BOTH_PVALUE_COL, BOTH_PVALUE_RANK, 
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    

TideLiteMatchSet::TideLiteMatchSet(ActivePeptideQueueLite* active_peptide_queue) {
  psm_scores_processed_ = false;
  active_peptide_queue_ = active_peptide_queue;
  psm_scores_ = PSMScores(active_peptide_queue->nPeptides_);  
};

TideLiteMatchSet::~TideLiteMatchSet() {
};

int* TideLiteMatchSet::getColumns(TSV_OUTPUT_FORMATS_T format, size_t& numHeaders){
  // TSV_OUTPUT_FORMATS_T is defined in ./src/model/objects.h
  switch (format) {
    case TIDE_SEARCH_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          numHeaders = sizeof(XCorr_tsv_cols) / sizeof(int);
          return XCorr_tsv_cols;
        case PVALUES:
          numHeaders = sizeof(Pvalues_tsv_cols) / sizeof(int);
          return Pvalues_tsv_cols;
        case DIAMETER:
          numHeaders = sizeof(Diameter_tsv_cols) / sizeof(int);
          return Diameter_tsv_cols;
      }
      break;
    case TIDE_SEARCH_MZTAB_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          numHeaders = sizeof(XCorr_mzTab_cols) / sizeof(int);
          return XCorr_mzTab_cols;
        // case PVALUES:
        //   numHeaders = sizeof(Pvalues_mzTab_cols) / sizeof(int);
        //   return Pvalues_mzTab_cols;
        // case DIAMETER:
        //   numHeaders = sizeof(Diameter_mzTab_cols) / sizeof(int);
        //   return Diameter_mzTab_cols;
      }
      break;
    case TIDE_SEARCH_PIN_TSV:
      break;
  }
}

string TideLiteMatchSet::getHeader(TSV_OUTPUT_FORMATS_T format, string tide_index_mztab_param_file) { 

  string header;  // Return variable
  size_t numHeaders = 0;
  int* header_cols = NULL;
  std::stringstream mztab_meta_data;
  std::string line;
  std::ifstream tide_index_mztab_param(tide_index_mztab_param_file);
  int param_cnt = 1;
  string mods;
  int cnt = 1;
  string site_prefix;
  string position_prefix;
  vector<string> tide_spectra_files;
  int search_engine_score_index = 1;    

  switch (format) {
    case TIDE_SEARCH_TSV:
      header_cols = getColumns(format, numHeaders);
      header += get_column_header(header_cols[0]);
      for (size_t i = 1; i < numHeaders; ++i) {
        header += '\t';
        header += get_column_header(header_cols[i]);
      }
      header += '\n';
      return header;
    case TIDE_SEARCH_MZTAB_TSV:

      mztab_meta_data << "MTD\tmzTab-version\t1.0.0\n";
      mztab_meta_data << "MTD\tmzTab-mode\tSummary\n";
      mztab_meta_data << "MTD\tmzTab-type\tIdentification\n";
      mztab_meta_data << "MTD\tdescription\tTide identification file.\n";
      tide_spectra_files = Params::GetStrings("tide spectra file");

      for(unsigned int i=0; i < tide_spectra_files.size(); i++) {
        mztab_meta_data << "MTD\tms_run[" << i + 1 << "]-location\t" << tide_spectra_files[i]  <<  "\n";
      }
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1001155, The SEQUEST result 'XCorr'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1001143, The SEQUEST result 'DeltaCn'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003358, XCorr rank]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003359, exact p-value'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003360, refactored XCorr'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003361, res-ev score'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003362, res-ev rank'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003363, res-ev p-value'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003364, combined p-value'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003365, combined p-value rank'.]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003366, tailor score'.]\n";
      
      // Add tide-index params
      if (tide_index_mztab_param) {
        while( std::getline(tide_index_mztab_param, line) ) {
          mztab_meta_data << line << "\n";
        }
      }
      // Add tide-search params
      mztab_meta_data << "MTD\tsoftware[2]\t[MS, MS:1002575, tide-search, "  << CRUX_VERSION << "]\n";    
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tmod-precision = " << Params::GetInt("mod-precision") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tauto-precursor-window = " << Params::GetString("auto-precursor-window") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tmax-precursor-charge = " << Params::GetInt("max-precursor-charge") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tprecursor-window = " << Params::GetDouble("precursor-window") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tprecursor-window-type = " << Params::GetString("precursor-window-type") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tauto-mz-bin-width = " << Params::GetString("auto-mz-bin-width") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tdeisotope = " << Params::GetDouble("deisotope") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tisotope-error = " << Params::GetString("isotope-error") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tmin-peaks = " << Params::GetInt("min-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tmz-bin-offset = " << Params::GetDouble("mz-bin-offset") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tmz-bin-width = " << Params::GetDouble("mz-bin-width") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tscore-function = " << Params::GetString("score-function") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tfragment-tolerance = " << Params::GetDouble("fragment-tolerance") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tevidence-granularity = " << Params::GetInt("evidence-granularity") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tremove-precursor-peak = " << (Params::GetBool("remove-precursor-peak")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tremove-precursor-tolerance = " << Params::GetDouble("remove-precursor-tolerance") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tscan-number = " << Params::GetString("scan-number") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tskip-preprocessing = " << (Params::GetBool("skip-preprocessing")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tspectrum-charge = " << Params::GetString("spectrum-charge") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tspectrum-max-mz = " << Params::GetDouble("spectrum-max-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tspectrum-min-mz = " << Params::GetDouble("spectrum-min-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tuse-flanking-peaks = " << (Params::GetBool("use-flanking-peaks")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tuse-neutral-loss-peaks = " << Params::GetDouble("use-neutral-loss-peaks") << "\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tnum-threads = " << Params::GetInt("num-threads") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-charges = " << Params::GetString("pm-charges") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-max-frag-mz = " << Params::GetDouble("pm-max-frag-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-max-precursor-delta-ppm = " << Params::GetDouble("pm-max-precursor-delta-ppm") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-max-precursor-mz = " << Params::GetDouble("pm-max-precursor-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-max-scan-separation = " << Params::GetInt("pm-max-scan-separation") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-min-common-frag-peaks = " << Params::GetInt("pm-min-common-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-min-frag-mz = " << Params::GetDouble("pm-min-frag-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-min-peak-pairs = " << Params::GetInt("pm-min-peak-pairs") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-min-precursor-mz = " << Params::GetDouble("pm-min-precursor-mz") <<"\n"; 
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-min-scan-frag-peaks = " << Params::GetInt("pm-min-scan-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-pair-top-n-frag-peaks = " << Params::GetInt("pm-pair-top-n-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tpm-top-n-frag-peaks = " << Params::GetInt("pm-top-n-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tconcat = " << (Params::GetBool("concat")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tspectrum-parser = " << Params::GetString("spectrum-parser") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\ttop-match = " << Params::GetInt("top-match") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-settings[" << param_cnt++ << "]\tuse-z-line = " << Params::GetInt("use-z-line") <<"\n";

      site_prefix = "";
      position_prefix = "Anywhere";
      mods = GetModificationList(MassConstants::mod_table_, site_prefix, position_prefix, false, cnt); //false = static mods;
      mztab_meta_data << mods;
      site_prefix = "N-term";
      position_prefix = "Any N-term";
      mods = GetModificationList(MassConstants::n_mod_table_, site_prefix, position_prefix, false, cnt); //false = static mods;
      mztab_meta_data << mods;
      site_prefix = "N-term";
      position_prefix = "Protein N-term";
      mods = GetModificationList(MassConstants::nprot_mod_table_, site_prefix, position_prefix, false, cnt); //false = static mods;
      mztab_meta_data << mods;
      site_prefix = "C-term";
      position_prefix = "Any C-term";
      mods = GetModificationList(MassConstants::c_mod_table_, site_prefix, position_prefix, false, cnt); //false = static mods;
      mztab_meta_data << mods;
      site_prefix = "C-term";
      position_prefix = "Protein C-term";
      mods = GetModificationList(MassConstants::cprot_mod_table_, site_prefix, position_prefix, false, cnt); //false = static mods;
      mztab_meta_data << mods;
      if ( cnt == 1 ) 
        mztab_meta_data << "MTD\tMS:1002453 (No fixed modifications searched)\n";

      cnt = 1;
      site_prefix = "";
      position_prefix = "Anywhere";
      mods = GetModificationList(MassConstants::mod_table_, site_prefix, position_prefix, true, cnt); //true = variable mods;
      mztab_meta_data << mods;
      site_prefix = "N-term";
      position_prefix = "Any N-term";
      mods = GetModificationList(MassConstants::n_mod_table_, site_prefix, position_prefix, true, cnt); //true = variable mods;
      mztab_meta_data << mods;
      site_prefix = "N-term";
      position_prefix = "Protein N-term";
      mods = GetModificationList(MassConstants::nprot_mod_table_, site_prefix, position_prefix, true, cnt); //true = variable mods;
      mztab_meta_data << mods;
      site_prefix = "C-term";
      position_prefix = "Any C-term";
      mods = GetModificationList(MassConstants::c_mod_table_, site_prefix, position_prefix, true, cnt); //true = variable mods;
      mztab_meta_data << mods;
      site_prefix = "C-term";
      position_prefix = "Protein C-term";
      mods = GetModificationList(MassConstants::cprot_mod_table_, site_prefix, position_prefix, true, cnt); //true = variable mods;
      mztab_meta_data << mods;
      if ( cnt == 1 ) 
        mztab_meta_data << "MTD\tMS:1002454 (No variable modifications searched)\m";
      header = mztab_meta_data.str();

      numHeaders = 0;
      header_cols = getColumns(format, numHeaders);
      header += get_column_header(header_cols[0]);
      for (size_t i = 1; i < numHeaders; ++i) {
        header += '\t';
        header += get_column_header(header_cols[i]);
      }
      header += '\n';
      return header;
    case TIDE_SEARCH_PIN_TSV:
      numHeaders = 0;
      header_cols = getColumns(format, numHeaders);
      header += get_column_header(header_cols[0]);
      for (size_t i = 1; i < numHeaders; ++i) {
        header += '\t';
        header += get_column_header(header_cols[i]);
      }
      header += '\n';
      return header;
  }

}

void TideLiteMatchSet::getReport(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, string &concat_or_target_report, string& decoy_report) { 

  // Clear the output variables  
  concat_or_target_report.clear();
  decoy_report.clear();

  // Get the top_n target and decoy PSMs
  gatherTargetsDecoys(); 
  
  //calculate tailor, delta_cn, and delta_lcn for the top n matches
  calculateAdditionalScores(concat_or_target_psm_scores_);  
  calculateAdditionalScores(decoy_psm_scores_);  // decoy_psm_scores is empty in case of concat=T

  // Prepare the results in a string
  printResults(format, spectrum_filename, sc, spectrum_file_cnt, concat_or_target_psm_scores_, concat_or_target_report);
  printResults(format, spectrum_filename, sc, spectrum_file_cnt, decoy_psm_scores_, decoy_report); // decoy_report is empty string if decoy_psm_scores is empty 

}

void TideLiteMatchSet::gatherTargetsDecoys() { 
  if (psm_scores_processed_ == true)
    return;

  // reset output variables
  concat_or_target_psm_scores_.clear();
  decoy_psm_scores_.clear();

  if (psm_scores_.size() == 0) {
    psm_scores_processed_ = true;
    return;
  }
  bool (*comp)(const Scores& x, const Scores& y);

  switch (curScoreFunction_) {
  case XCORR_SCORE:
    comp = &cmpXcorrScore;
    break;
  case PVALUES:
    comp = &cmpCombinedPvalue;
    break;
  case DIAMETER:
    break;
  }

  quantile_score_ = 1.0;

  // Calculate Tailor scores. Get the 99th quantile:
  // if (curScoreFunction_ == XCORR_SCORE) {

  int quantile_pos = (int)(TAILOR_QUANTILE_TH*(double)psm_scores_.size()+0.5)-1; // zero indexed
  if (quantile_pos < 2) 
    quantile_pos = 2;  // the third element
  if (quantile_pos >= psm_scores_.size()) 
    quantile_pos = psm_scores_.size()-1; // the last element

  make_heap(psm_scores_.begin(), psm_scores_.end(), cmpXcorrScore);
  for (int i = 0; i <= quantile_pos; ++i){
    pop_heap(psm_scores_.begin(), psm_scores_.end()-i, cmpXcorrScore);
    Scores back = psm_scores_[psm_scores_.size()-1-i];
  }
  quantile_score_ = psm_scores_[psm_scores_.size()-1-quantile_pos].xcorr_score_ +TAILOR_OFFSET; // Make sure scores positive
  // }

  // get the value of the last score for the delta_lcn scores
  last_psm_ = std::min_element(psm_scores_.begin(), psm_scores_.end(), comp);
  // printf("last value:%lf\n", (*last_psm_).xcorr_score_);
  
  // Gather target and decoy PSMs
  int gatherSize = top_matches_ + 1; // Get one more psms than the top-matches, so the delta_cn can be calculated correctly for the last element.
  map<int, int> decoyWriteCount;
  
  make_heap(psm_scores_.begin(), psm_scores_.end(), comp);

  for (PSMScores::iterator i = psm_scores_.end(); i > psm_scores_.begin(); ) {
    pop_heap(psm_scores_.begin(), i--, comp);
    if ((*i).active_ == false)
      continue;

// TODO: uncommented for debugging, put it back
    PeptideLite* peptide = active_peptide_queue_->GetPeptide((*i).ordinal_);
    if (concat_ || !peptide->IsDecoy()) {
      if (concat_or_target_psm_scores_.size() < gatherSize) {
        concat_or_target_psm_scores_.push_back(*i);
      }
    } else {
      int idx = peptide->DecoyIdx();
      map<int, int>::iterator j = decoyWriteCount.find(idx);
      if (j == decoyWriteCount.end()) {
        j = decoyWriteCount.insert(make_pair(idx, 0)).first;
      }
      if (j->second < gatherSize) {
        j->second++;
        decoy_psm_scores_.push_back(*i);
      }
    }

    if ((concat_ && concat_or_target_psm_scores_.size() >= gatherSize) || (!concat_ && concat_or_target_psm_scores_.size() >= gatherSize && decoy_psm_scores_.size() >= gatherSize*decoy_num_)){
      break;
    }
  }  
}

void TideLiteMatchSet::calculateAdditionalScores(PSMScores& psm_scores) {  // Additional scores are:  delta_cn, delta_lcn, tailor;
  // The  gatherTargetsDecoys must be run before calling this function.

  switch (curScoreFunction_) {
  case XCORR_SCORE:
  case DIAMETER:
    last_psm_ = psm_scores.end()-1;
    for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it){
      (*it).tailor_ = ((*it).xcorr_score_  + TAILOR_OFFSET )/ quantile_score_;
      (*it).delta_lcn_ = ((*it).xcorr_score_ - (*last_psm_).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = ((*it).xcorr_score_ - (*(it+1)).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      else 
        (*it).delta_cn_ = 0.0;
      // printf("xcorr:%lf, tailor:%lf, delta_cn: %lf, delta_lcn: %lf\n", (*it).xcorr_score_, (*it).tailor_, (*it).delta_cn_, (*it).delta_lcn_);
    }
    break;
  case PVALUES:
    last_psm_ = psm_scores.end()-1;
    for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it){
      (*it).tailor_ = ((*it).xcorr_score_  + TAILOR_OFFSET )/ quantile_score_;
      (*it).delta_lcn_ = -log10((*it).combined_pval_) + log10((*last_psm_).combined_pval_);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = -log10((*it).combined_pval_) + log10((*(it+1)).combined_pval_);
      else 
        (*it).delta_cn_ = 0.0;
      // printf("xcorr:%lf, tailor:%lf, delta_cn: %lf, delta_lcn: %lf, ", (*it).xcorr_score_, (*it).tailor_, (*it).delta_cn_, (*it).delta_lcn_);
      // printf("refactored xcorr: %lf, p-value: %lf, CPV: %lf\n ", (*it).refactored_xcorr_, (*it).exact_pval_, (*it).combinedPval_, (*it).delta_lcn_);
    }
    break;
  case PVALUES_HR:
  case PVALUES_LR:
  case HYPERSCORE:
  case HYPERSCORE_LA:

    break;
  }

  psm_scores_processed_ = true;

}

void TideLiteMatchSet::printResults(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, PSMScores& psm_scores, string& report) { 
  // The order of the fields of the results is solely based on the column order
  size_t numHeaders;
  int* header_cols = getColumns(format, numHeaders);
  int cnt = 1;
  for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it, ++cnt) {
    if (cnt > top_matches_)
      break;
    PeptideLite* peptide = active_peptide_queue_->GetPeptide((*it).ordinal_);

    string proteinNames = peptide->GetLocationStr(decoy_prefix_);
    string flankingAAs = peptide->GetFlankingAAs();
    string peptide_with_mods = peptide->SeqWithMods(mod_precision_);
    string crux_modifications;
    string mztab_modifications;
    peptide->getModifications(mod_precision_, crux_modifications, mztab_modifications); 
    string peptide_without_mods = peptide->Seq();   
    // The order of the fields depends on the columns defined at the beginning of this file.
    // The fields below can be in arbitrary order.
    for (size_t i = 0; i < numHeaders; ++i) {   
      switch (header_cols[i]){
      case FILE_COL:
        report += spectrum_filename;
        break;
      case SCAN_COL:
        report += StringUtils::ToString(sc->spectrum->SpectrumNumber(), 0); // Scan Id
        break;
      case MZTAB_CHARGE:
      case CHARGE_COL:
        report += StringUtils::ToString(sc->charge, 0);                     // Charge
        break;
      case MZTAB_EXP_MASS_TO_CHARGE:
      case SPECTRUM_PRECURSOR_MZ_COL:
        report += StringUtils::ToString(sc->spectrum->PrecursorMZ(), mass_precision_);                             //spectrum precursor mz
        break;
      case MZTAB_OPT_MS_RUN_1_SPECTRUM_NEUTRAL_MASS:
      case SPECTRUM_NEUTRAL_MASS_COL:
        report += StringUtils::ToString((sc->spectrum->PrecursorMZ() - MASS_PROTON)*sc->charge, mass_precision_);  // spectrum neutral mass
        break;
      case MZTAB_CALC_MASS_TO_CHARGE:
      case PEPTIDE_MASS_COL:
        report += StringUtils::ToString(peptide->Mass(), mass_precision_);                                         // spectrum neutral mass
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_2:
      case DELTA_CN_COL:
        report += StringUtils::ToString((*it).delta_cn_, score_precision_);         // delta_cn
        break;
      case MZTAB_OPT_MS_RUN_1_DELTA_LCN:
      case DELTA_LCN_COL:
        report += StringUtils::ToString((*it).delta_lcn_, score_precision_);        // delta_lcn
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_1:
      case XCORR_SCORE_COL:
        report += StringUtils::ToString((*it).xcorr_score_, score_precision_);      // xcorr score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_5:
      case REFACTORED_SCORE_COL:
        report += StringUtils::ToString((*it).refactored_xcorr_, score_precision_);      // refactored xcorr score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_4:
      case EXACT_PVALUE_COL:
        report += StringUtils::ToString((*it).exact_pval_, score_precision_, false);      // exact p-value score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_6:
      case RESIDUE_EVIDENCE_COL:
        report += StringUtils::ToString((*it).resEv_score_, score_precision_);      // res-ev score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_8:
      case RESIDUE_PVALUE_COL:
        report += StringUtils::ToString((*it).resEv_pval_, score_precision_, false);      // res_ev pvalue
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_9:
      case BOTH_PVALUE_COL:
        report += StringUtils::ToString((*it).combined_pval_, score_precision_, false);      // combined p-value score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_11:
      case TAILOR_COL:
        report += StringUtils::ToString((*it).tailor_, score_precision_);           // tailor score
        break;
      case BY_IONS_MATCHED_COL:
        report += StringUtils::ToString((*it).by_ion_matched_, score_precision_);   // by ions matched
        break;
      case BY_IONS_TOTAL_COL:
        report += StringUtils::ToString((*it).by_ion_total_, score_precision_);     // total no of by ions
        break;
      case BY_IONS_FRACTION_COL:
        report += StringUtils::ToString((double)((*it).by_ion_matched_)/(*it).by_ion_total_, score_precision_);  // fraction of the matched per total by-ions
        break;
      case BY_IONS_REPEAT_MATCH_COL:
        report += StringUtils::ToString((*it).repeat_ion_match_, score_precision_);  // fraction of the matched per total by-ions
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_3:
      case MZTAB_SEARCH_ENGINE_SCORE_7:
      case MZTAB_SEARCH_ENGINE_SCORE_10:
      case BOTH_PVALUE_RANK:
      case XCORR_RANK_COL:
        report += StringUtils::ToString(cnt, 0);
        break;
      case MZTAB_OPT_MS_RUN_1_DISTINCT_MATCHES_PER_SPEC:
      case DISTINCT_MATCHES_SPECTRUM_COL:
        report += StringUtils::ToString(active_peptide_queue_->nCandPeptides_, 0); //StringUtils::ToString(target ? n_concat_or_target_matches_:n_decoy_matches_, 0);  // no candidate peptides
        break;
      case SEQUENCE_COL:
        report += peptide_with_mods;        // peptide sequence with modifications
        break;
      case MODIFICATIONS_COL:
        report += crux_modifications;            // the list of modifications in the peptide
        break;
      case MZTAB_SEQUENCE:                  // MZTAB column
      case UNMOD_SEQUENCE_COL:
        report += peptide_without_mods;     // plain peptide sequence stripped with mods
        break;
      case MZTAB_ACCESSION:                 // MZTAB column
      case PROTEIN_ID_COL:
        report += proteinNames;             // protein IDs
        break;
      case FLANKING_AA_COL:
        report += flankingAAs;              // flanking Amino Acids 
        break;
      case MZTAB_OPT_MS_RUN_1_TARGET_DECOY:
      case TARGET_DECOY_COL:
        report += string(peptide->IsDecoy()? "decoy":"target");  // target or decoy
        break;
      case MZTAB_OPT_MS_RUN_1_ORIGINAL_TARGET_SEQUENCE_COL:
      case ORIGINAL_TARGET_SEQUENCE_COL:
        report += peptide->TargetSeq();     // original target sequences
        break;
      case MZTAB_OPT_MS_RUN_1_DECOY_INDEX:
      case DECOY_INDEX_COL:
        report += StringUtils::ToString(peptide->DecoyIdx(), 0);  // decoy index
        break;
 // MZTAB columns      
      case MZTAB_PSH:
        report += "PSM";  // Just prints PSM
        break;
      case MZTAB_PSM_ID:
        report += StringUtils::ToString(psm_id_mzTab_++, 0);  // PSM ID
        break;
      case MZTAB_UNIQUE:
      case MZTAB_DATABASE_VERSION:
        report += "null";  // null
        break;
      case MZTAB_DATABASE:
        report += fasta_file_name_;  // the fasta file name
        break;
      case MZTAB_SEARCH_ENGINE:
        report += "software[2]";  // the fasta file name
        break;
      case MZTAB_MODIFICATIONS:
        report += mztab_modifications;  // the fasta file name
        break;
      case RETENTION_TIME_COL:
      case MZTAB_RETENTION_TIME:
        report += StringUtils::ToString(sc->spectrum->RTime(), mass_precision_);  ;  // retention time
        break;
      case MZTAB_SPECTRA_REF:
        report += "ms_run["+ StringUtils::ToString(spectrum_file_cnt, 0) +"]:index=" +  StringUtils::ToString(sc->spectrum->SpectrumNumber(), 0);  // the fasta file name
        break;
      case MZTAB_PRE:
        report += flankingAAs.substr(0, 1);  // the fasta file name
        break;
      case MZTAB_POST:
        report += flankingAAs.substr(1, 1);  // the fasta file name
        break;
      case MZTAB_START:
        report += StringUtils::ToString(peptide->FirstLocPos(), 0);  // the fasta file name
        break;
      case MZTAB_END:
        report += StringUtils::ToString(peptide->FirstLocPos()+peptide->Len()-1, 0);  // the fasta file name
        break;
      }

      if (i < numHeaders-1)  // If not the last column, add a column separator
        report += '\t';
    }
    report += '\n'; 
  }
}

string TideLiteMatchSet::GetModificationList(const pb::ModTable* mod_table, string site_prefix, string position_prefix, bool variable_mods, int& cnt){
  string mods;
  string site;
  if (variable_mods == false) { 
    for (int i = 0; i < mod_table->static_mod_size(); i++) {
      
      const pb::Modification& mod = mod_table->static_mod(i);

      if (mod.has_name() && mod.has_amino_acids()) {
        site = mod.amino_acids();
        if (site == "X")
          site = "";
        mods += "MTD\tfixed_mod[" + std::to_string(cnt) + "]\t["   + mod.name()            + "]\n";
        mods += "MTD\tfixed_mod[" + std::to_string(cnt) + "]-site\t" + site_prefix + site    + "\n";
        mods += "MTD\tfixed_mod[" + std::to_string(cnt++)+"]-position\t" + position_prefix + "\n";
      }
    }
  } else {
    for (int i = 0; i < mod_table->variable_mod_size(); i++) {
      
      const pb::Modification& mod = mod_table->variable_mod(i);

      if (mod.has_name() && mod.has_amino_acids()) {
        site = mod.amino_acids();
        if (site == "X")
          site = "";
        mods += "MTD\tvariable_mod[" + std::to_string(cnt) + "]\t["     +  mod.name()           + "]\n";
        mods += "MTD\tvariable_mod[" + std::to_string(cnt) + "]-site\t" + site_prefix + site    + "\n";
        mods += "MTD\tvariable_mod[" + std::to_string(cnt++)+"]-position\t" + position_prefix + "\n";
      }
    }

  }
  return mods;
}
