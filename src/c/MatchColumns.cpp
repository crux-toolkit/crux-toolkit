/*************************************************************************//**
 * \file MatchColumns.cpp
 * \brief Just keeps track of column names for match files.
 ****************************************************************************/

#include "MatchColumns.h"
#include "carp.h"

static const char* match_column_strings[NUMBER_MATCH_COLUMNS] = {
  "scan",
  "charge",
  "spectrum precursor m/z",
  "spectrum neutral mass",
  "peptide mass",
  "delta_cn",
  "sp score",
  "sp rank",
  "xcorr score",
  "xcorr rank",
  "p-value",
#ifdef NEW_COLUMNS
  "Weibull PSM q-value",
  "Weibull peptide q-value",    // NEW
  "decoy PSM q-value",
  "decoy peptide q-value",      // NEW
  "percolator score",
  "percolator rank",
  "percolator PSM q-value",
  "percolator peptide q-value", // NEW
  "q-ranker score",
  "q-ranker PSM q-value",
  "q-ranker peptide q-value",   // NEW
#else
  "Weibull est. q-value",
  "decoy q-value (xcorr)",
  "percolator score",
  "percolator rank",
  "percolator q-value",
  "q-ranker score",
  "q-ranker q-value",
#endif
  "b/y ions matched",
  "b/y ions total",
  "matches/spectrum",
  "sequence",
  "cleavage type",
  "protein id",
  "flanking aa",
  "unshuffled sequence",
  "eta",
  "beta",
  "shift",
  "corr",
  "SIN",
  "NSAF",
  "EMPAI",
  "parsimony rank"
};

/**
 * Get the name of a given column, by index.
 */
const char* get_column_header(
  int columnIndex
) {
  if ((columnIndex < 0) || (columnIndex >= NUMBER_MATCH_COLUMNS)) {
    carp(CARP_FATAL, "Cannot access output column %d.\n", columnIndex);
  }
  return(match_column_strings[columnIndex]);
}
