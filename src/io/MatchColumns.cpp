/*************************************************************************
 * \file MatchColumns.cpp
 * \brief Just keeps track of column names for match files.
 *************************************************************************/

/**************************************************************************
 * Steps for adding a new column to Crux tab-delimited output files.
 *
 * o Add the header string for the column to the match_column_strings
 *   list in MatchColumns.cpp.
 *
 * o Add a corresponding const for the column index to enum
 *   MATCH_COLUMNS_T in MatchColumns.h
 *
 * o Add logic to addColumnNames in io/MatchFileWriter.cpp to ensure
 *   that the column gets included for the appropriate commands.
 *
 * o If necessary, add the column to the headers array in
 *   app/TideMatchSet.cpp and the cols_to_print array in
 *   app/AssignConfidenceApplication.cpp.
 *
 * o If necessary, modify the writeToFile function in
 *   app/TideMatchSet.cpp to print the new column.  Be sure to modify
 *   both variants of the function.
 *
 * o Update the smoke tests to include the new column in the stored
 *   results files.
 *
 * o Update doc/user/file-formats/txt-format.html to be sure that
 *   every column listed in MatchColumns.h appears there, in the
 *   correct order.
 **************************************************************************/

#include "MatchColumns.h"
#include "io/carp.h"
#include <string.h>

static const char* match_column_strings[NUMBER_MATCH_COLUMNS] = {
  "file",
  "file_idx",
  "scan",
  "charge",
  "retention time",  
  "spectrum precursor m/z",
  "spectrum neutral mass",
  "peptide mass",
  "delta_cn",
  "delta_lcn",
  "sp score",
  "sp rank",
  "xcorr score",
  "xcorr rank",
  "exact p-value",
  "refactored xcorr",
  "res-ev p-value", //Added by Andy Lin
  "res-ev score", //Added by Andy Lin
  "res-ev rank", //Added by Andy Lin
  "combined p-value", //Added by Andy Lin
  "combined p-value rank", //Added by Andy Lin
  "Sidak adjusted p-value",  
  "tailor score",  //Added for tailor score calibration method by AKF  
  "e-value",
  "smoothed p-value",
  "distinct matches/peptide",
  "precursor intensity logrank M0", // added by Yang
  "precursor intensity logrank M1", // added by Yang
  "precursor intensity logrank M2", // added by Yang
  "rt-diff", // added by Yang
  "dynamic fragment p-value", // added by Yang
  "static fragment p-value", // added by Yang
  "precursor coelution", // added by Yang
  "fragment coelution", // added by Yang
  "precursor fragment coelution", // added by Yang
  "ensemble score", // added by Yang
#ifdef NEW_COLUMNS
  "decoy PSM q-value",
  "decoy peptide q-value",      // NEW
  "percolator score",
  "percolator rank",
  "percolator PSM q-value",
  "percolator peptide q-value", // NEW
#else
  "decoy q-value (xcorr)",
  "decoy PEP (xcorr)",
  "decoy q-value (e-value)",
  "decoy PEP (e-value)",
  "percolator score",
  "percolator rank",
  "percolator q-value",
  "percolator PEP",
#endif
  "tdc q-value",
  "mix-max q-value",
  "b/y ions matched",
  "b/y ions total",
  "b/y ions fraction",
  "b/y ion repeat match",
  "total matches/spectrum",
  "distinct matches/spectrum",
  "sequence",
  "modifications",
  "cleavage type",
  "unmodified sequence",  
  "protein id",
  "peptides",
  "flanking aa",
  "target/decoy",
  "original target sequence",
  "RAW",
  "SIN",
  "NSAF",
  "dNSAF",
  "EMPAI",
  "parsimony rank",
  "decoy matches/spectrum",
  "SpecId",  // for PinWriter, many of these are repeats with slightly different names, is there better way to do this?
  "Label",
  "ScanNr",
  "ExpMass",
  "CalcMass",
  "lnrSp",
  "deltLCn",
  "deltCn",
  "Xcorr",
  "Sp",
  "IonFrac",
  "Mass",
  "PepLen",
  "enzN",
  "enzC",
  "enzInt",
  "lnNumSP",
  "dm",
  "absdM",
  "Peptide",
  "Proteins", // end for PinWriter.
  "ppm",
  "xcorr 1",
  "xcorr 2",
  "protein id x",
  "index name",
  "decoy index",
  // Percolator PIN columns.
  "PSMId",
  "score",
  "q-value",
  "posterior_error_prob",
  "peptide",
  "proteinIds",
 // mzTab Headers
  "PSH",
  "sequence",
  "PSM_ID",
  "accession",
  "unique",
  "database",
  "database_version",
  "search_engine",
  "search_engine_score[1]",   // [MS, MS:1001155, The SEQUEST result 'XCorr']
  "search_engine_score[2]",  // [MS, MS:1003366, tailor score]
  "search_engine_score[3]",   // [MS, MS:1001143, The SEQUEST result 'DeltaCn']
  "search_engine_score[4]",   // [MS, MS:1003358, XCorr rank]
  "search_engine_score[5]",   // [MS, MS:1003360, refactored XCorr]
  "search_engine_score[6]",   // [MS, MS:1003359, exact p-value]
  "search_engine_score[7]",   // [MS, MS:1003361, res-ev score]
  "search_engine_score[8]",   // [MS, MS:1003363, res-ev p-value]
  "search_engine_score[9]",   // [MS, MS:1003364, combined p-value]
  "search_engine_score[10]",  // [MS, MS:1003365, combined p-value rank]
  "search_engine_score[11]",  // [MS, MS:1002354, PSM-level q-value]
  "modifications",
  "retention_time",
  "charge",
  "exp_mass_to_charge",
  "calc_mass_to_charge",
  "spectra_ref",
  "pre",
  "post",
  "start",
  "end",
  "opt_ms_run[1]_spectrum_neutral_mass",
  "opt_ms_run[1]_delta_lcn",
  "opt_ms_run[1]_distinct_matches_per_spectrum",
  "opt_ms_run[1]_target_or_decoy",
  "opt_ms_run[1]_original_target_sequence",
  "opt_ms_run[1]_decoy_index",

};

/**
 * Get the name of a given column, by index.
 */
const char* get_column_header(
  int columnIndex
) {
  if (columnIndex < 0 || columnIndex >= NUMBER_MATCH_COLUMNS) {
    carp(CARP_FATAL, "Cannot access output column %d.\n", columnIndex);
  }
  carp(CARP_DETAILED_DEBUG, "get_column_header: %d/%d %s",
       columnIndex, NUMBER_MATCH_COLUMNS - 1, match_column_strings[columnIndex]);
  return(match_column_strings[columnIndex]);
}

int get_column_idx(
  const char* column_name
) {
  for (int i = FILE_COL; i < NUMBER_MATCH_COLUMNS; ++i) {
    if (strcmp(column_name, match_column_strings[i]) == 0) {
      return i;
    }
  }
  return INVALID_COL;
}

