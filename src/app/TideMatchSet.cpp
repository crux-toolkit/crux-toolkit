#include <iomanip>

#include "TideMatchSet.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "tide/peptide.h"
#include "crux_version.h"
#include "TideSearchApplication.h"

// SCORE_FUNCTION_T is defined in ./src/model/objects.h
SCORE_FUNCTION_T TideMatchSet::curScoreFunction_ = INVALID_SCORE_FUNCTION;
int TideMatchSet::top_matches_ = 0;
int TideMatchSet::decoy_num_ = 0;
int TideMatchSet::mass_precision_ = 0;
int TideMatchSet::score_precision_ = 0;
int TideMatchSet::mod_precision_ = 0;
bool TideMatchSet::concat_ = false;
string TideMatchSet::decoy_prefix_ = "";
int TideMatchSet::psm_id_mzTab_  = 1;
string TideMatchSet::fasta_file_name_ = "null";


// column IDs are defined in ./src/io/MatchColumns.h and /src/io/MatchColumns.cpp
// The order of the columns in the output is defined here by the order of the column Ids
int TideMatchSet::XCorr_tsv_cols[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, RETENTION_TIME_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
    XCORR_RANK_COL, DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
 int TideMatchSet::Pvalues_tsv_cols[] = {  //TODO: update the columns.
    FILE_COL, SCAN_COL, CHARGE_COL, RETENTION_TIME_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
    RESIDUE_EVIDENCE_COL, RESIDUE_PVALUE_COL, BOTH_PVALUE_COL, BOTH_PVALUE_RANK, 
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    

  int TideMatchSet::Diameter_tsv_cols[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL,  XCORR_SCORE_COL, TAILOR_COL, XCORR_RANK_COL,
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
    PRECURSOR_INTENSITY_RANK_M0_COL, PRECURSOR_INTENSITY_RANK_M1_COL, PRECURSOR_INTENSITY_RANK_M2_COL,
    RT_DIFF_COL, DYN_FRAGMENT_PVALUE_COL, STA_FRAGMENT_PVALUE_COL,
    COELUTE_MS1_COL, COELUTE_MS2_COL, COELUTE_MS1_MS2_COL, ENSEMBLE_SCORE_COL,
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL
  };

int TideMatchSet::XCorr_mzTab_cols[] = {
    MZTAB_PSH, MZTAB_SEQUENCE, MZTAB_PSM_ID, MZTAB_ACCESSION, MZTAB_UNIQUE, MZTAB_DATABASE,
    MZTAB_DATABASE_VERSION, MZTAB_SEARCH_ENGINE, 
    MZTAB_SEARCH_ENGINE_SCORE_1,  // [MS, MS:1001155, The SEQUEST result 'XCorr'.]
    MZTAB_SEARCH_ENGINE_SCORE_2,  // [MS, MS:1003366, tailor score'.]
    MZTAB_SEARCH_ENGINE_SCORE_3,  // [MS, MS:1001143, The SEQUEST result 'DeltaCn'.]
    MZTAB_SEARCH_ENGINE_SCORE_4,  // [MS, MS:1003358, XCorr rank]
    // MZTAB_SEARCH_ENGINE_SCORE_5,  // [MS, MS:1003360, refactored XCorr'.]
    // MZTAB_SEARCH_ENGINE_SCORE_6,  // [MS, MS:1003359, exact p-value'.]
    // MZTAB_SEARCH_ENGINE_SCORE_7,  // [MS, MS:1003361, res-ev score'.]
    // MZTAB_SEARCH_ENGINE_SCORE_8,  // [MS, MS:1003363, res-ev p-value'.]
    // MZTAB_SEARCH_ENGINE_SCORE_9,  // [MS, MS:1003364, combined p-value'.]
    // MZTAB_SEARCH_ENGINE_SCORE_10, // [MS, MS:1003365, combined p-value rank'.]
    MZTAB_MODIFICATIONS, MZTAB_RETENTION_TIME,
    MZTAB_CHARGE, MZTAB_EXP_MASS_TO_CHARGE, MZTAB_CALC_MASS_TO_CHARGE, MZTAB_SPECTRA_REF,
    MZTAB_PRE, MZTAB_POST, MZTAB_START, MZTAB_END, MZTAB_OPT_MS_RUN_1_SPECTRUM_NEUTRAL_MASS,
    MZTAB_OPT_MS_RUN_1_DELTA_LCN, MZTAB_OPT_MS_RUN_1_DISTINCT_MATCHES_PER_SPEC,
    MZTAB_OPT_MS_RUN_1_TARGET_DECOY, MZTAB_OPT_MS_RUN_1_ORIGINAL_TARGET_SEQUENCE_COL, 
    MZTAB_OPT_MS_RUN_1_DECOY_INDEX
  };    

int TideMatchSet::Pvalues_mzTab_cols[] = {
    MZTAB_PSH, MZTAB_SEQUENCE, MZTAB_PSM_ID, MZTAB_ACCESSION, MZTAB_UNIQUE, MZTAB_DATABASE,
    MZTAB_DATABASE_VERSION, MZTAB_SEARCH_ENGINE, 
    MZTAB_SEARCH_ENGINE_SCORE_1,  // [MS, MS:1001155, The SEQUEST result 'XCorr'.]
    MZTAB_SEARCH_ENGINE_SCORE_2,  // [MS, MS:1003366, tailor score'.]
    MZTAB_SEARCH_ENGINE_SCORE_3,  // [MS, MS:1001143, The SEQUEST result 'DeltaCn'.]
    // MZTAB_SEARCH_ENGINE_SCORE_4,  // [MS, MS:1003358, XCorr rank]
    MZTAB_SEARCH_ENGINE_SCORE_5,  // [MS, MS:1003360, refactored XCorr'.]
    MZTAB_SEARCH_ENGINE_SCORE_6,  // [MS, MS:1003359, exact p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_7,  // [MS, MS:1003361, res-ev score'.]
    MZTAB_SEARCH_ENGINE_SCORE_8,  // [MS, MS:1003363, res-ev p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_9,  // [MS, MS:1003364, combined p-value'.]
    MZTAB_SEARCH_ENGINE_SCORE_10, // [MS, MS:1003365, combined p-value rank'.]
    MZTAB_MODIFICATIONS, MZTAB_RETENTION_TIME,
    MZTAB_CHARGE, MZTAB_EXP_MASS_TO_CHARGE, MZTAB_CALC_MASS_TO_CHARGE, MZTAB_SPECTRA_REF,
    MZTAB_PRE, MZTAB_POST, MZTAB_START, MZTAB_END, MZTAB_OPT_MS_RUN_1_SPECTRUM_NEUTRAL_MASS,
    MZTAB_OPT_MS_RUN_1_DELTA_LCN, MZTAB_OPT_MS_RUN_1_DISTINCT_MATCHES_PER_SPEC,
    MZTAB_OPT_MS_RUN_1_TARGET_DECOY, MZTAB_OPT_MS_RUN_1_ORIGINAL_TARGET_SEQUENCE_COL, 
    MZTAB_OPT_MS_RUN_1_DECOY_INDEX
  };    

// int TideMatchSet::XCorr_pin_cols[] = {
//     POUT_PSMID_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
//     PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
//     BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
//     XCORR_RANK_COL, DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
//     PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
//     DECOY_INDEX_COL
//   };    
//  int TideMatchSet::Pvalues_pin_cols[] = {  //TODO: update the colums.
//     FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
//     PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
//     RESIDUE_EVIDENCE_COL, RESIDUE_PVALUE_COL, BOTH_PVALUE_COL, BOTH_PVALUE_RANK, 
//     DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
//     PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
//     DECOY_INDEX_COL
//   };    

TideMatchSet::TideMatchSet(ActivePeptideQueue* active_peptide_queue, ObservedPeakSet* observed) {
  psm_scores_processed_ = false;
  active_peptide_queue_ = active_peptide_queue;
  observed_ = observed;  // Pointer to the experimental spectrum data 

  psm_scores_ = PSMScores(active_peptide_queue->nPeptides_);  
};

TideMatchSet::~TideMatchSet() {
};

int* TideMatchSet::getColumns(TSV_OUTPUT_FORMATS_T format, size_t& numHeaders){
  // TSV_OUTPUT_FORMATS_T is defined in ./src/model/objects.h
  switch (format) {
    case DIAMETER_TSV:
      numHeaders = sizeof(Diameter_tsv_cols) / sizeof(int);
      return Diameter_tsv_cols;
    case TIDE_SEARCH_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          numHeaders = sizeof(XCorr_tsv_cols) / sizeof(int);
          return XCorr_tsv_cols;
        case PVALUES:
          numHeaders = sizeof(Pvalues_tsv_cols) / sizeof(int);
          return Pvalues_tsv_cols;
        // case DIAMETER:
        //   numHeaders = sizeof(Diameter_tsv_cols) / sizeof(int);
        //   return Diameter_tsv_cols;
      }
      break;
    case TIDE_SEARCH_MZTAB_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          numHeaders = sizeof(XCorr_mzTab_cols) / sizeof(int);
          return XCorr_mzTab_cols;
        case PVALUES:
          numHeaders = sizeof(Pvalues_mzTab_cols) / sizeof(int);
          return Pvalues_mzTab_cols;
        // case DIAMETER:
        //   numHeaders = sizeof(Diameter_mzTab_cols) / sizeof(int);
        //   return Diameter_mzTab_cols;
      }
      break;
    case TIDE_SEARCH_PIN_TSV:  // Consider to print pin format directly, so that MakePinApplication can be removed.
      break;
  }
}

string TideMatchSet::getHeader(TSV_OUTPUT_FORMATS_T format, string tide_index_mztab_param_file) { 

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
    case DIAMETER_TSV:
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
        mztab_meta_data << "MTD\tms_run[" << i + 1 << "]-location\tfile://" << tide_spectra_files[i]  <<  "\n";
      }
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1001155, SEQUEST:xcorr, ]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003366, tailor score, ]\n";
      mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1001143, SEQUEST:deltacn, ]\n";
      if (curScoreFunction_ != PVALUES) {
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003358, XCorr rank, ]\n";
      }
      if (curScoreFunction_ == PVALUES) {
        search_engine_score_index++; //a jamp for the XCorr rank
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003360, refactored XCorr, ]\n";
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003359, exact p-value, ]\n";
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003361, res-ev score, ]\n";
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003363, res-ev p-value, ]\n";
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003364, combined p-value, ]\n";
        mztab_meta_data << "MTD\tpsm_search_engine_score[" << search_engine_score_index++ << "]\t[MS, MS:1003365, combined p-value rank, ]\n";
      }
      // Add tide-index params
      if (tide_index_mztab_param) {
        while( std::getline(tide_index_mztab_param, line) ) {
          mztab_meta_data << line << "\n";
        }
      }
      // Add tide-search params
      mztab_meta_data << "MTD\tsoftware[2]\t[MS, MS:1002575, tide-search, "  << CRUX_VERSION << "]\n";    
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tmod-precision = " << Params::GetInt("mod-precision") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tauto-precursor-window = " << Params::GetString("auto-precursor-window") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tmax-precursor-charge = " << Params::GetInt("max-precursor-charge") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tprecursor-window = " << Params::GetDouble("precursor-window") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tprecursor-window-type = " << Params::GetString("precursor-window-type") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tauto-mz-bin-width = " << Params::GetString("auto-mz-bin-width") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tdeisotope = " << Params::GetDouble("deisotope") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tisotope-error = " << Params::GetString("isotope-error") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tmin-peaks = " << Params::GetInt("min-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tmz-bin-offset = " << Params::GetDouble("mz-bin-offset") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tmz-bin-width = " << Params::GetDouble("mz-bin-width") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tscore-function = " << Params::GetString("score-function") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tfragment-tolerance = " << Params::GetDouble("fragment-tolerance") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tevidence-granularity = " << Params::GetInt("evidence-granularity") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tremove-precursor-peak = " << (Params::GetBool("remove-precursor-peak")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tremove-precursor-tolerance = " << Params::GetDouble("remove-precursor-tolerance") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tscan-number = " << Params::GetString("scan-number") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tskip-preprocessing = " << (Params::GetBool("skip-preprocessing")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tspectrum-charge = " << Params::GetString("spectrum-charge") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tspectrum-max-mz = " << Params::GetDouble("spectrum-max-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tspectrum-min-mz = " << Params::GetDouble("spectrum-min-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tuse-flanking-peaks = " << (Params::GetBool("use-flanking-peaks")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tuse-neutral-loss-peaks = " << Params::GetDouble("use-neutral-loss-peaks") << "\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tnum-threads = " << Params::GetInt("num-threads") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-charges = " << Params::GetString("pm-charges") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-max-frag-mz = " << Params::GetDouble("pm-max-frag-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-max-precursor-delta-ppm = " << Params::GetDouble("pm-max-precursor-delta-ppm") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-max-precursor-mz = " << Params::GetDouble("pm-max-precursor-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-max-scan-separation = " << Params::GetInt("pm-max-scan-separation") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-min-common-frag-peaks = " << Params::GetInt("pm-min-common-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-min-frag-mz = " << Params::GetDouble("pm-min-frag-mz") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-min-peak-pairs = " << Params::GetInt("pm-min-peak-pairs") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-min-precursor-mz = " << Params::GetDouble("pm-min-precursor-mz") <<"\n"; 
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-min-scan-frag-peaks = " << Params::GetInt("pm-min-scan-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-pair-top-n-frag-peaks = " << Params::GetInt("pm-pair-top-n-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tpm-top-n-frag-peaks = " << Params::GetInt("pm-top-n-frag-peaks") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tconcat = " << (Params::GetBool("concat")?"True":"False") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tspectrum-parser = " << Params::GetString("spectrum-parser") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\ttop-match = " << Params::GetInt("top-match") <<"\n";
      mztab_meta_data << "MTD\tsoftware[2]-setting[" << param_cnt++ << "]\tuse-z-line = " << Params::GetInt("use-z-line") <<"\n";

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
        mztab_meta_data << "MTD\tMS:1002454 (No variable modifications searched)\n";
      header = mztab_meta_data.str();

      numHeaders = 0;
      header_cols = getColumns(format, numHeaders);  // numHeaders is an output variable
      header += get_column_header(header_cols[0]);
      for (size_t i = 1; i < numHeaders; ++i) {
        header += '\t';
        header += get_column_header(header_cols[i]);
      }
      header += '\n';
      return header;
    case TIDE_SEARCH_PIN_TSV:
    // TODO: Finish Percolator Input (PIN) File format later
/*      numHeaders = 0;
      header_cols = getColumns(format, numHeaders);
      header += get_column_header(header_cols[0]);
      for (size_t i = 1; i < numHeaders; ++i) {
        header += '\t';
        header += get_column_header(header_cols[i]);
      }
      header += '\n';
*/
      return header;
  }
}

void TideMatchSet::getReport(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, string &concat_or_target_report, string& decoy_report) { 

  // Clear the output variables  
  concat_or_target_report.clear();
  decoy_report.clear();

  // Get the top_n target and decoy PSMs
  gatherTargetsDecoys(); 
  
  //calculate tailor, delta_cn, and delta_lcn for the top n matches
  calculateAdditionalScores(concat_or_target_psm_scores_, sc);  
  calculateAdditionalScores(decoy_psm_scores_, sc);  // decoy_psm_scores is empty in case of concat=T

  // Prepare the results in a string
  printResults(format, spectrum_filename, sc, spectrum_file_cnt, true, concat_or_target_psm_scores_, concat_or_target_report);  // true = target
  printResults(format, spectrum_filename, sc, spectrum_file_cnt, false, decoy_psm_scores_, decoy_report); // decoy_report is an empty string if decoy_psm_scores is empty; false = decoy

}

void TideMatchSet::gatherTargetsDecoys() { 
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
  // case DIAMETER:
  case XCORR_SCORE:
    comp = &cmpXcorrScore;
    break;
  case PVALUES:
    comp = &cmpCombinedPvalue;
    break;
  }

  quantile_score_ = 1.0;

  // Calculate Tailor scores. Get the 99th quantile:

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
  quantile_score_ = psm_scores_[psm_scores_.size()-1-quantile_pos].xcorr_score_ + TAILOR_OFFSET; // Make sure scores positive

  // get the value of the last score for the delta_lcn scores
  last_psm_ = std::min_element(psm_scores_.begin(), psm_scores_.end(), comp);
  
  // Gather target and decoy PSMs
  int gatherSize = top_matches_ + 1; // Get one more psms than the top-matches, so the delta_cn can be calculated correctly for the last rankedd PSM element.
  map<int, int> decoyWriteCount;
  
  make_heap(psm_scores_.begin(), psm_scores_.end(), comp);

  for (PSMScores::iterator i = psm_scores_.end(); i > psm_scores_.begin(); ) {
    pop_heap(psm_scores_.begin(), i--, comp);
    if ((*i).active_ == false)
      continue;

    Peptide* peptide = active_peptide_queue_->GetPeptide((*i).ordinal_);
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

void TideMatchSet::calculateAdditionalScores(PSMScores& psm_scores, const SpectrumCollection::SpecCharge* sc) {  // Additional scores are:  delta_cn, delta_lcn, tailor;
  // The  gatherTargetsDecoys must be run before calling this function.
  int last_psm_pos = -2;
  if (top_matches_ >= psm_scores.size()) {
    last_psm_pos = -1;
  }
  last_psm_ = psm_scores.end()+last_psm_pos;
  
  int temp;
  int repeat_ion_match;  
  Peptide* peptide;

  for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it){
    // Count the repeating matching ions. This was used in SP scoring
    temp = 0;
    repeat_ion_match = 0;
    peptide = (*((*it).peptide_itr_));
    temp = TideSearchApplication::PeakMatching(*observed_, peptide->peaks_1b, temp, repeat_ion_match);
    temp = TideSearchApplication::PeakMatching(*observed_, peptide->peaks_1y, temp, repeat_ion_match);

    if (sc->charge > 2) {
      temp = TideSearchApplication::PeakMatching(*observed_, peptide->peaks_2b, temp, repeat_ion_match);
      temp = TideSearchApplication::PeakMatching(*observed_, peptide->peaks_2y, temp, repeat_ion_match);
    }
    (*it).repeat_ion_match_ = repeat_ion_match;

    // Perform Tailor calibration
    (*it).tailor_ = ((*it).xcorr_score_  + TAILOR_OFFSET )/ quantile_score_;

    switch (curScoreFunction_) {
    case XCORR_SCORE:
      (*it).delta_lcn_ = ((*it).xcorr_score_ - (*last_psm_).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = ((*it).xcorr_score_ - (*(it+1)).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      else 
        (*it).delta_cn_ = 0.0;
      break;
    case PVALUES:
      (*it).delta_lcn_ = -log10((*it).combined_pval_) + log10((*last_psm_).combined_pval_);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = -log10((*it).combined_pval_) + log10((*(it+1)).combined_pval_);
      else 
        (*it).delta_cn_ = 0.0;
    }
    // }
    // break;
  // case PVALUES_HR:
  // case PVALUES_LR:
  // case HYPERSCORE:
  // case HYPERSCORE_LA:
  }

  psm_scores_processed_ = true;

}

void TideMatchSet::printResults(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, const SpectrumCollection::SpecCharge* sc, int spectrum_file_cnt, bool target, PSMScores& psm_scores, string& report,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* intensity_map,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* logrank_map,
    map<PSMScores::iterator, boost::tuple<double, double, double>>* coelute_map,
    map<PSMScores::iterator, boost::tuple<double, double>>* ms2pval_map,
    map<string, double>* peptide_predrt_map) { 
  // The order of the fields of the results is solely based on the column order
  size_t numHeaders;
  int* header_cols = getColumns(format, numHeaders);

/*
The following counter called cnt, is used to count the PSMs reported. The decoy index of a
target peptide is -1, and a decoy index of a decoy peptide used to be 0. Later, somebody else
added an option to make tide-search handle multiple decoy sets (This was developed by Uri for
the Average TDC method). So, the decoy index has become a zero-based index indicating which 
decoy set this given decoy peptide is originating from.

Tide-search  supports concatenated search and separated search and it is also able to report 
the top-N matches per spectrum. Multiple decoy sets are supported only in separated 
target-decoy search. 

In the case of only target or concatenated target decoy search, the decoy index can be 0 or 1,
and we use the cnt[0] to count how many target and decoys PSMs have been reported so far. 

In the case of separated target decoy search and this function reports only targets, then cnt[0] 
counts only targets.
In the case of separated target decoy search and this function reports only decoys, then the 
cnt[i] counts only decoys, for i = 0-->decoy_num 
*/

  vector<int> cnt;
  for (int i=0; i <= decoy_num_; ++i) {
    cnt.push_back(0);
  }
  double predrt;
  
  for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it) {
    Peptide* peptide = active_peptide_queue_->GetPeptide((*it).ordinal_);
    int decoy_idx = peptide->DecoyIdx();
    decoy_idx = decoy_idx < 0 ? 0 : decoy_idx;
    ++cnt[decoy_idx];
    if (cnt[decoy_idx] > top_matches_) {
      continue;
    }

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
      case MZTAB_SEARCH_ENGINE_SCORE_3:
      case DELTA_CN_COL:
        report += StringUtils::ToString((*it).delta_cn_, score_precision_);         // delta_cn
        break;
      case MZTAB_OPT_MS_RUN_1_DELTA_LCN:
      case DELTA_LCN_COL:
        report += StringUtils::ToString((*it).delta_lcn_, score_precision_);        // delta_lcn
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_1:  // xcorr score
      case XCORR_SCORE_COL:
        report += StringUtils::ToString((*it).xcorr_score_, score_precision_);      // xcorr score
        break;       
      case MZTAB_SEARCH_ENGINE_SCORE_5:   // refactored XCorr
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString((*it).refactored_xcorr_, score_precision_, false);      // refactored XCorr
        } else {
          report += "null";       // refactored XCorr
        }
        break;
      case REFACTORED_SCORE_COL:
        report += StringUtils::ToString((*it).refactored_xcorr_, score_precision_);      // refactored xcorr score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_6:      // exact p-value
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString((*it).exact_pval_, score_precision_, false);      // exact p-value score
        } else {
          report += "null";      // exact p-value score
        }
        break;
      case EXACT_PVALUE_COL:
        report += StringUtils::ToString((*it).exact_pval_, score_precision_, false);      // exact p-value score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_7:                                                 // res-ev score
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString((*it).resEv_score_, score_precision_, false);      // res-ev score
        } else {
          report += "null";      // res-ev score
        }
        break;
      case RESIDUE_EVIDENCE_COL:
        report += StringUtils::ToString((*it).resEv_score_, score_precision_);      // res-ev score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_8:     // res-ev p-value rank
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString((*it).resEv_pval_, score_precision_, false);      // res-ev score
        } else {
          report += "null";      // exact p-value score
        }
        break;
      case RESIDUE_PVALUE_COL:
        report += StringUtils::ToString((*it).resEv_pval_, score_precision_, false);      // res_ev pvalue
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_9:  //  combined p-value'
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString((*it).combined_pval_, score_precision_, false);      // combined p-value
        } else {
          report += "null";      // combined p-value
        }
        break;
      case BOTH_PVALUE_COL:
        report += StringUtils::ToString((*it).combined_pval_, score_precision_, false);      // combined p-value score
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_2:
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
  
      case MZTAB_SEARCH_ENGINE_SCORE_4:  // [MS, MS:1003358, XCorr rank]
        if (curScoreFunction_ != PVALUES) {
          report += StringUtils::ToString(cnt[decoy_idx], 0);  // rank
        } else {
          report += "null";      // xcorr rank in p-value scoring
        }
        break;
      case MZTAB_SEARCH_ENGINE_SCORE_10:  // [MS, MS:1003365, combined p-value rank'.]
        if (curScoreFunction_ == PVALUES) {
          report += StringUtils::ToString(cnt[decoy_idx], 0);  // rank
        } else {
          report += "null";      // exact p-value score
        }
        break;
      case BOTH_PVALUE_RANK:    // combined p-value rank
      case XCORR_RANK_COL:
        report += StringUtils::ToString(cnt[decoy_idx], 0);
        break;
      case MZTAB_OPT_MS_RUN_1_DISTINCT_MATCHES_PER_SPEC:
      case DISTINCT_MATCHES_SPECTRUM_COL:
        if (concat_ == true) {
          report += StringUtils::ToString(active_peptide_queue_->nCandPeptides_, 0); // Print num of targets and decoys
        } else if (target) {
          report += StringUtils::ToString(active_peptide_queue_->CandPeptidesTarget_, 0); // Print num of targets
        } else {
          report += StringUtils::ToString(active_peptide_queue_->CandPeptidesDecoy_, 0); // Print num of decoys
        }
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
        report += "PSM";  // Just prints PSM at the beginnning of the row
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
        report += "[MS, MS:1002575, tide-search, ]";  // the fasta file name
        break;
      case MZTAB_MODIFICATIONS:
        report += mztab_modifications;  // the fasta file name
        break;
      case RETENTION_TIME_COL:
      case MZTAB_RETENTION_TIME:
        report += StringUtils::ToString(sc->spectrum->RTime(), mass_precision_);  ;  // retention time
        break;
      case MZTAB_SPECTRA_REF:
        report += "ms_run["+ StringUtils::ToString(spectrum_file_cnt+1, 0) +"]:index=" +  StringUtils::ToString(sc->spectrum->SpectrumNumber(), 0);  // the fasta file name
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
      // Diameter features: PRECURSOR_INTENSITY_RANK_M0_COL, PRECURSOR_INTENSITY_RANK_M1_COL, PRECURSOR_INTENSITY_RANK_M2_COL,
      // RT_DIFF_COL, DYN_FRAGMENT_PVALUE_COL, STA_FRAGMENT_PVALUE_COL,
      // COELUTE_MS1_COL, COELUTE_MS2_COL, COELUTE_MS1_MS2_COL, ENSEMBLE_SCORE_COL,
      case PRECURSOR_INTENSITY_RANK_M0_COL:
        if (intensity_map != NULL) {
          boost::tuple<double, double, double> intensity_tuple = intensity_map->at(it);
          report += StringUtils::ToString(intensity_tuple.get<0>()+intensity_tuple.get<1>()+intensity_tuple.get<2>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case PRECURSOR_INTENSITY_RANK_M1_COL:
        if (intensity_map != NULL) {
          boost::tuple<double, double, double> intensity_tuple = intensity_map->at(it);
          report += StringUtils::ToString(intensity_tuple.get<0>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case PRECURSOR_INTENSITY_RANK_M2_COL:
        if (logrank_map != NULL) {
          boost::tuple<double, double, double> logrank_tuple = logrank_map->at(it);
          report += StringUtils::ToString(logrank_tuple.get<0>()+logrank_tuple.get<1>()+logrank_tuple.get<2>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case RT_DIFF_COL:
        predrt = 0.5;
        if (peptide_predrt_map != NULL) {
          map<string, double>::iterator predrtIter = peptide_predrt_map->find(peptide_with_mods);
          if (predrtIter != peptide_predrt_map->end()) { 
            predrt = predrtIter->second; 
          }
        }
        report += StringUtils::ToString(fabs(predrt - sc->spectrum->RTime()), score_precision_, true);
        break;
      case DYN_FRAGMENT_PVALUE_COL:
        if (ms2pval_map != NULL) {
          boost::tuple<double, double> ms2pval = ms2pval_map->at(it);
          report += StringUtils::ToString(ms2pval.get<0>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case STA_FRAGMENT_PVALUE_COL:
        if (ms2pval_map != NULL) {
          boost::tuple<double, double> ms2pval = ms2pval_map->at(it);
          report += StringUtils::ToString(ms2pval.get<1>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case COELUTE_MS1_COL:
        if (coelute_map != NULL) {
          boost::tuple<double, double, double> coelute_tuple = coelute_map->at(it);
          report += StringUtils::ToString(coelute_tuple.get<0>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case COELUTE_MS2_COL:
        if (coelute_map != NULL) {
          boost::tuple<double, double, double> coelute_tuple = coelute_map->at(it);
          report += StringUtils::ToString(coelute_tuple.get<1>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case COELUTE_MS1_MS2_COL:
        if (coelute_map != NULL) {
          boost::tuple<double, double, double> coelute_tuple = coelute_map->at(it);
          report += StringUtils::ToString(coelute_tuple.get<2>(), score_precision_, true);
        } else {
          report += "0";
        }
        break;
      case ENSEMBLE_SCORE_COL:
        report += StringUtils::ToString(0.0, score_precision_, true);
        break;
      }

      if (i < numHeaders-1)  // If not the last column, add a column separator
        report += '\t';
    }
    report += '\n'; 
  }
}

/*
  int TideMatchSet::Diameter_tsv_cols[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL,  XCORR_SCORE_COL, TAILOR_COL, XCORR_RANK_COL,
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, BY_IONS_REPEAT_MATCH_COL,
    PRECURSOR_INTENSITY_RANK_M0_COL, PRECURSOR_INTENSITY_RANK_M1_COL, PRECURSOR_INTENSITY_RANK_M2_COL,
    RT_DIFF_COL, DYN_FRAGMENT_PVALUE_COL, STA_FRAGMENT_PVALUE_COL,
    COELUTE_MS1_COL, COELUTE_MS2_COL, COELUTE_MS1_MS2_COL, ENSEMBLE_SCORE_COL,
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL
  };
*/
string TideMatchSet::GetModificationList(const pb::ModTable* mod_table, string site_prefix, string position_prefix, bool variable_mods, int& cnt){
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
