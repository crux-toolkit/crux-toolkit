/**
 * \file MatchFileWriter.cpp
 * DATE: October 26, 2010
 * AUTHOR: Barbara Frewen
 * \brief Object for writing tab-deliminted text files of PSMs (matches).
 * This class is the extension of DelimitedFileWriter where the
 * columns are those in MatchColumns.  Which columns are written to
 * file depend on the COMMAND_TYPE_T.  For some commands, they also
 * depend on the columns found in the input file as given by a
 * MatchFileReader. 
 */

#include "MatchFileWriter.h"
#include "parameter.h"
#include "util/Params.h"
#include "app/TideSearchApplication.h"
#include <iostream>

using namespace std;

/**
 * \returns A blank MatchFileWriter object.
 */
MatchFileWriter::MatchFileWriter() 
  : DelimitedFileWriter(),
    num_columns_(0) {
  for(int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    match_to_print_[col_type] = false;
    match_precision_[col_type] = 0;
    match_fixed_float_[col_type] = true;
  }
  setPrecision();
} 

/**
 * \returns A blank MatchFileWriter object and opens a file for
 * writing.
 */
MatchFileWriter::MatchFileWriter(const char* filename) 
  : DelimitedFileWriter(filename),
    num_columns_(0) {
  for(int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    match_to_print_[col_type] = false;
    match_precision_[col_type] = 0;
    match_fixed_float_[col_type] = true;
  }
  setPrecision();
}

/**
 * Destructor
 */
MatchFileWriter::~MatchFileWriter() {
}

/**
 * Set the correct level of precision for each MATCH_COLUMNS_T type.
 * Also set whether the field should be fixed float or not.
 */
void MatchFileWriter::setPrecision() {
  for(int col_idx = 0; col_idx < NUMBER_MATCH_COLUMNS; col_idx++) {
    switch(col_idx) {
      // integer and string fields
    case FILE_COL:
    case SCAN_COL:
    case CHARGE_COL:
    case SP_RANK_COL:
    case XCORR_RANK_COL:
    case PERCOLATOR_RANK_COL:
    case BY_IONS_MATCHED_COL:
    case BY_IONS_TOTAL_COL:
    case DISTINCT_MATCHES_SPECTRUM_COL:
    case MATCHES_SPECTRUM_COL:
    case SEQUENCE_COL:
    case CLEAVAGE_TYPE_COL:
    case PROTEIN_ID_COL:
    case FLANKING_AA_COL:
    case ORIGINAL_TARGET_SEQUENCE_COL:
    case PARSIMONY_RANK_COL:
    case RAW_SCORE_COL:  //Raw counts should be integral
      match_precision_[col_idx] = 0;
      match_fixed_float_[col_idx] = true;
      break;

      // mass fields
    case SPECTRUM_PRECURSOR_MZ_COL:
    case SPECTRUM_NEUTRAL_MASS_COL:
    case MASS_COL:
    case DM_COL:
    case ABS_DM_COL:
    case PEPTIDE_MASS_COL:
      match_precision_[col_idx] = Params::GetInt("mass-precision");
      match_fixed_float_[col_idx] = true;
      break;

      // score fields
    case DELTA_CN_COL:
    case SP_SCORE_COL:
    case XCORR_SCORE_COL:
    case XCORR_FIRST_COL:
    case XCORR_SECOND_COL:
    case EVALUE_COL:
    case PVALUE_COL:
    case WEIBULL_QVALUE_COL:
    case WEIBULL_PEP_COL:
    case DECOY_XCORR_QVALUE_COL:
    case DECOY_XCORR_PEP_COL:
    case DECOY_EVALUE_QVALUE_COL:
    case DECOY_EVALUE_PEP_COL:
    case PERCOLATOR_SCORE_COL:
    case PERCOLATOR_QVALUE_COL:
    case PERCOLATOR_PEP_COL:
    case QRANKER_SCORE_COL:
    case QRANKER_QVALUE_COL:
    case QRANKER_PEP_COL:
    case SIN_SCORE_COL:
    case NSAF_SCORE_COL:
    case DNSAF_SCORE_COL:
    case EMPAI_SCORE_COL:
    case LNR_SP_COL:
    case DELT_L_CN_COL:
    case DELT_CN_COL:
    case XCORR_COL:
    case SP_COL:
    case ION_FRAC_COL:
    case LN_NUM_SP_COL:
    case EXACT_PVALUE_COL:
    case REFACTORED_SCORE_COL:
    case RESIDUE_EVIDENCE_COL:
    case BOTH_PVALUE_COL:
    case SIDAK_ADJUSTED_COL:
    case QVALUE_MIXMAX_COL:
    case QVALUE_TDC_COL:
#ifdef NEW_COLUMNS
    case WEIBULL_PEPTIDE_QVALUE_COL:      // NEW
    case DECOY_XCORR_PEPTIDE_QVALUE_COL:  // NEW
    case PERCOLATOR_PEPTIDE_QVALUE_COL:   // NEW
    case QRANKER_PEPTIDE_QVALUE_COL:      // NEW
#endif
      match_precision_[col_idx] = Params::GetInt("precision");
      match_fixed_float_[col_idx] = false;
      break;

      // special cases
    case ETA_COL:
    case BETA_COL:
    case SHIFT_COL:
    case CORR_COL:
      match_precision_[col_idx] = 6;
      match_fixed_float_[col_idx] = false;
      break;

    case NUMBER_MATCH_COLUMNS:
    case INVALID_COL:
      carp(CARP_FATAL, "Invalid match column type for setting precision.");
      break;
    }
  }
}

/**
 * Defines the columns to print based on the vector of flags
 * indiciating if the MATCH_COLUMN_T should be printed.
 */
void MatchFileWriter::addColumnNames(const std::vector<bool>& col_is_printed) {
  // for each column, if we should print it, mark as true
  for(size_t col_idx = 0; col_idx < col_is_printed.size(); col_idx++) {
    bool print_it = col_is_printed[col_idx];
    if( print_it ) {
      match_to_print_[col_idx] = true;
    } 
  }
}

/**
 * Adds another columns to print.  Printed in order the names are set.
 */
void MatchFileWriter::addColumnName(MATCH_COLUMNS_T column_type) {
  match_to_print_[column_type] = true;
}

/**
 * Adds which columns to print based on the COMMAND_TYPE_T. Only for
 * search-for-matches, sequest-search and spectral-counts.
 */
void MatchFileWriter::addColumnNames(CruxApplication* application, 
                                     bool has_decoys) {

  COMMAND_T command = application->getCommand();

  switch (command) {
  // commands with no tab files
  case MISC_COMMAND:
  case PROCESS_SPEC_COMMAND: ///< print-processed-spectra
  case GENERATE_PEPTIDES_COMMAND: ///< generate-peptides
  case GET_MS2_SPECTRUM_COMMAND: ///< get-ms2-spectrum 
  case PREDICT_PEPTIDE_IONS_COMMAND: ///< predict-peptide-ions
  case VERSION_COMMAND:      ///< just print the version number
  // invalid
  case NUMBER_COMMAND_TYPES:
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command (%s) for creating a MatchFileWriter.",
         application->getName().c_str());
    return;

  // commands that also require list of cols to print
  case QVALUE_COMMAND:       ///< compute-q-values
  case PERCOLATOR_COMMAND:
  case QRANKER_COMMAND:
  case BARISTA_COMMAND:
    carp(CARP_FATAL, 
         "Post-search command %s requires a list of columns to print.",
         application->getName().c_str());
    return;

  // valid commands
  case TIDE_SEARCH_COMMAND: ///< tide-search
    if (Params::GetBool("file-column")) {
      addColumnName(FILE_COL);
    }
    if (Params::GetBool("compute-sp") || Params::GetBool("sqt-output")) {
      if (Params::GetBool("exact-p-value")) {
        addColumnName(EXACT_PVALUE_COL);
        addColumnName(REFACTORED_SCORE_COL);
      } else {
        addColumnName(SP_SCORE_COL);
      }
      addColumnName(SP_RANK_COL);
      addColumnName(BY_IONS_MATCHED_COL);
      addColumnName(BY_IONS_TOTAL_COL);
    }
    break;

  case LOCALIZE_MODIFICATION_COMMAND:
    if (Params::GetBool("file-column")) {
      addColumnName(FILE_COL);
    }
    addColumnName(SCAN_COL);
    addColumnName(CHARGE_COL);
    addColumnName(SPECTRUM_PRECURSOR_MZ_COL);
    addColumnName(SPECTRUM_NEUTRAL_MASS_COL);
    addColumnName(PEPTIDE_MASS_COL);
    addColumnName(XCORR_SCORE_COL);
    addColumnName(SEQUENCE_COL);
    addColumnName(MODIFICATIONS_COL);
    addColumnName(PROTEIN_ID_COL);
    addColumnName(FLANKING_AA_COL);
    addColumnName(TARGET_DECOY_COL);
    return;

  case XLINK_SEARCH_COMMAND:
    if (Params::GetBool("compute-p-values")) {
      addColumnName(PVALUE_COL);
      addColumnName(ETA_COL);
      addColumnName(BETA_COL);
      addColumnName(SHIFT_COL);
      addColumnName(CORR_COL);
    }
    if (!Params::GetBool("use-old-xlink")) {
      if (Params::GetInt("xlink-top-n") != 0) {
        addColumnName(XCORR_FIRST_COL);
        addColumnName(XCORR_SECOND_COL);
      }
      if (Params::GetBool("file-column")) {
        addColumnName(FILE_COL);
      }
      addColumnName(ENZ_INT_COL);
    }
    addColumnName(XLINK_TYPE_COL);
    break;

  case SPECTRAL_COUNTS_COMMAND:
    // protein or peptide
    if( string_to_quant_level_type(Params::GetString("quant-level")) == PEPTIDE_QUANT_LEVEL ) {
      addColumnName(SEQUENCE_COL);
    } else {
      addColumnName(PROTEIN_ID_COL);
      // parsimony?
      if( string_to_parsimony_type(Params::GetString("parsimony")) != PARSIMONY_NONE ) {
        addColumnName(PARSIMONY_RANK_COL);
      }
    }

    // SIN or NSAF score
    switch (string_to_measure_type(Params::GetString("measure"))) {
      case MEASURE_RAW:
        addColumnName(RAW_SCORE_COL);
        break;
      case MEASURE_SIN:
        addColumnName(SIN_SCORE_COL);
        break;
      case MEASURE_NSAF:
        addColumnName(NSAF_SCORE_COL);
        break;
      case MEASURE_DNSAF:
        addColumnName(DNSAF_SCORE_COL);
        break;
      case MEASURE_EMPAI:
        addColumnName(EMPAI_SCORE_COL);
        break;
      default:
        break; //do nothing  
    }
    return; // do not add additional columns
  }

  // All search commands have these columns
  addColumnName(SCAN_COL);
  addColumnName(CHARGE_COL);
  addColumnName(SPECTRUM_PRECURSOR_MZ_COL);
  addColumnName(SPECTRUM_NEUTRAL_MASS_COL);
  addColumnName(PEPTIDE_MASS_COL);
  addColumnName(DELTA_CN_COL);
  if (Params::GetBool("exact-p-value")) {
    addColumnName(EXACT_PVALUE_COL);
    addColumnName(REFACTORED_SCORE_COL);
  } else {
    addColumnName(XCORR_SCORE_COL);
  }
  addColumnName(XCORR_RANK_COL);
  addColumnName(DISTINCT_MATCHES_SPECTRUM_COL);
  addColumnName(SEQUENCE_COL);
  addColumnName(MODIFICATIONS_COL);
  addColumnName(CLEAVAGE_TYPE_COL);
  addColumnName(PROTEIN_ID_COL);
  addColumnName(FLANKING_AA_COL);
  addColumnName(TARGET_DECOY_COL);
  if ((has_decoys || Params::GetBool("concat")) &&
      !TideSearchApplication::proteinLevelDecoys()) {
    addColumnName(ORIGINAL_TARGET_SEQUENCE_COL);
  }

}

/**
 * Adds which columns to print based on the COMMAND_TYPE_T and a list
 * of columns to print. For all post-search commands.
 */
void MatchFileWriter::addColumnNames
  (CruxApplication* application, 
  bool has_decoys,
  const vector<bool>& cols_to_print) {

  COMMAND_T command = application->getCommand();

  switch (command) {
  // commands with no tab files
  case MISC_COMMAND:
  case PROCESS_SPEC_COMMAND: ///< print-processed-spectra
  case GENERATE_PEPTIDES_COMMAND: ///< generate-peptides
  case GET_MS2_SPECTRUM_COMMAND: ///< get-ms2-spectrum 
  case PREDICT_PEPTIDE_IONS_COMMAND: /// predict-peptide-ions
  case VERSION_COMMAND:      ///< just print the version number
  // invalid
  case SPECTRAL_COUNTS_COMMAND:
  case NUMBER_COMMAND_TYPES:
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command (%s) for creating a MatchFileWriter.",
         application->getName().c_str());
    return;

  // search commands handled elsewhere
  case XLINK_SEARCH_COMMAND:
    addColumnNames(application, has_decoys);
    return;

  // valid commands
  case QVALUE_COMMAND:       ///< compute-q-values
    if( cols_to_print[PVALUE_COL] ) {
      addColumnName(WEIBULL_QVALUE_COL);
//      addColumnName(WEIBULL_PEP_COL);
      //addColumnName(WEIBULL_PEPTIDE_QVALUE_COL);
    } else {
      if ( cols_to_print[EVALUE_COL]) {
//        addColumnName(DECOY_EVALUE_QVALUE_COL);
//        addColumnName(DECOY_EVALUE_PEP_COL);
      }
//     addColumnName(DECOY_XCORR_QVALUE_COL);
//     addColumnName(DECOY_XCORR_PEP_COL);
      //addColumnName(DECOY_XCORR_PEPTIDE_QVALUE_COL);
    }
    break;

  case PERCOLATOR_COMMAND:
    addColumnName(PERCOLATOR_SCORE_COL);
    addColumnName(PERCOLATOR_RANK_COL);
    addColumnName(PERCOLATOR_QVALUE_COL);
    addColumnName(PERCOLATOR_PEP_COL);
    break;

  case QRANKER_COMMAND:
    addColumnName(QRANKER_SCORE_COL);
    addColumnName(QRANKER_QVALUE_COL);
    addColumnName(QRANKER_PEP_COL);
    break;

  case BARISTA_COMMAND:
    // place holder
    break;
  }

  if( has_decoys ) {
    addColumnName(ORIGINAL_TARGET_SEQUENCE_COL);
  }

  // now add remaining columns from the input file
  addColumnNames(cols_to_print);

  // FIXME (BF 10-27-10): where do these go?
  //  PERCOLATOR_PEPTIDE_QVALUE_COL,   
  //  QRANKER_PEPTIDE_QVALUE_COL,      

}

/**
 * Write header to file using column names that have been set.
 */
void MatchFileWriter::writeHeader() {
  num_columns_ = 0;
  // set file position index for all columns being printed
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    if( match_to_print_[col_type] == true ) {
      match_indices_[col_type] = num_columns_++;
    } else {
      match_indices_[col_type] = -1;
    }
  }

  // set all the names for which we have match_indices_
  column_names_.assign(num_columns_, "");
  for(unsigned int col_type = 0; col_type < NUMBER_MATCH_COLUMNS; col_type++) {
    if( match_indices_[col_type] > -1 ) {
      if (get_column_header(col_type) == NULL) {
        carp(CARP_FATAL, "Error col type: %d doesn't exist!", col_type);
      }
      DelimitedFileWriter::setColumnName(get_column_header(col_type), 
                                         match_indices_[col_type]);
    }
  }

  DelimitedFileWriter::writeHeader();

  // every line will be this length, reserve space in current row for speed
  current_row_.assign(num_columns_, "");
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
