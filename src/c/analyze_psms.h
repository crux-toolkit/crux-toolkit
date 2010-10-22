/**
 * \file analyze_psms.h
 * The starting point for all crux post-search operations.
 */
#ifndef ANALYZE_PSMS_H
#define ANALYZE_PSMS_H

#include "objects.h"
#include "carp.h"
#include "OutputFiles.h"
#include "q-value.h"
#include "PercolatorCInterface.h"
#include "QRankerCInterface.h"

/**
 * \brief Takes a directory containing PSM files and a protein index
 * and analyzes the PSMs using compute-q-values, percolator or q-ranker.
 */
void analyze_matches_main(
  COMMAND_T command,
  int argc,
  char** argv
);

/**
 * \brief Analyze matches using the percolator or qranker algorithm.
 * 
 * Runs the specified algorithm on the PSMs in the psm_result_folder
 * for a search against the sequence database fasta_file. Optionally 
 * puts the percolator PSM feature vectors into feature_file, if it is 
 * not NULL.
 * \returns a pointer to a MATCH_COLLECTION_T object
 * \callgraph
 */
MATCH_COLLECTION_T* run_percolator_or_qranker(
  COMMAND_T command,                                          
  char* input_directory, 
  char* fasta_file, 
  OutputFiles& output);

#endif //ANALYZE_PSMS_H

