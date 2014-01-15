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
 * \returns a pointer to a MatchCollection object
 * \callgraph
 */
MatchCollection* run_percolator(
  char* input_directory, 
  char* fasta_file, 
  OutputFiles& output);

/**
 * Compute posterior error probabilities (PEP) from the given target
 * and decoy scores.
 * \returns A newly allocated array of PEP for the target scores
 * sorted.
 */
double* compute_PEP(double* target_scores, ///< scores for target matches
                    int num_targets,       ///< size of target_scores
                    double* decoy_scores,  ///< scores for decoy matches
                    int num_decoys,       ///< size of decoy_scores
                    bool ascending = false); ///< are the scores ascending/descending?

#endif //ANALYZE_PSMS_H

