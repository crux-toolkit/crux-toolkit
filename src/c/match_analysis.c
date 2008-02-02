/*****************************************************************************
 * \file match_analysis.c
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an 
 *              optional parameter file, search all the spectrum against 
 *              the peptides in the sequence database, and return the high 
 *              scoring peptides. 
 * REVISION: $ $
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "protein.h"
#include "peptide.h"
#include "spectrum.h"
#include "parse_arguments.h" //delete
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "scorer.h"
#include "match.h"
#include "match_collection.h"
#include "PercolatorCInterface.h"

#define MAX_PSMS 10000000
// 14th decimal place
#define EPSILON 0.00000000000001 
#define NUM_ANALYSIS_OPTIONS 5
#define NUM_ANALYSIS_ARGUMENTS 2

/** 
 * Routines to run various match analyses. Explained in more detail below.
 */
MATCH_COLLECTION_T* run_percolator(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file); 

MATCH_COLLECTION_T* run_qvalue(
  char* psm_result_folder, 
  char* fasta_file 
  ); 

MATCH_COLLECTION_T* run_nothing(
  char* psm_result_folder, 
  char* fasta_file,
  char* feature_file 
  ); 

int output_matches(
    MATCH_COLLECTION_T* match_collection,
    SCORER_TYPE_T scorer_type
    );

/***********************************************************************/
int main(int argc, char** argv){

  /* Declarations */
  char* psm_result_folder = NULL;
  char* fasta_file = NULL; //rename
  char* feature_file = NULL;

  /* Define command line arguments */
  int num_options = NUM_ANALYSIS_OPTIONS;
  char* option_list[NUM_ANALYSIS_OPTIONS] = {
    "verbosity",
    "parameter-file",
    "algorithm",
    "feature-file",
    "use-index" //not yet implemented, below set to true
  };

  int num_arguments = NUM_ANALYSIS_ARGUMENTS;
  char* argument_list[NUM_ANALYSIS_ARGUMENTS] = {
    "psm-folder",
    "protein input",
  };

  /* for debugging handling of parameters*/
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  set_verbosity_level(CARP_ERROR);

  /* Set up parameters and set defaults in parameter.c */
  initialize_parameters();

  /* Define optional and required arguments in parameter.c */
  select_cmd_line_options(option_list, num_options );
  select_cmd_line_arguments(argument_list, num_arguments);

  /* Parse the command line and optional paramter file
     does sytnax, type, and bounds checking and dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux-analyze-matches");

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /* Get arguments */
  psm_result_folder = get_string_parameter("psm-folder");
  fasta_file = get_string_parameter("protein input");
  feature_file = get_string_parameter("feature-file");

  /* Get options */
  ALGORITHM_TYPE_T algorithm_type = get_algorithm_type_parameter("algorithm");
  SCORER_TYPE_T scorer_type = PERCOLATOR_SCORE;
  MATCH_COLLECTION_T* match_collection = NULL;

  /* Perform the analysis */
  switch(algorithm_type){
  case PERCOLATOR_ALGORITHM:
    carp(CARP_INFO, "Running percolator");
    match_collection = run_percolator(psm_result_folder,
				      fasta_file,
				      feature_file);
    scorer_type = PERCOLATOR_SCORE;
    break;
    
  case QVALUE_ALGORITHM:
    carp(CARP_INFO, "Running qvalue");
    match_collection = run_qvalue(psm_result_folder, fasta_file);
    scorer_type = Q_VALUE;
    break;
    
  case NO_ALGORITHM:
    carp(CARP_INFO, "No analysis algorithm chosen.");
    match_collection = run_nothing(psm_result_folder,
				   fasta_file,
				   feature_file);
    scorer_type = XCORR; // TODO put in something to default to the primary
    // score in the run
    break;
  default:
    ;
  }  
  
  carp(CARP_INFO, "Outputting matches.");
  output_matches(match_collection, scorer_type);
  // MEMLEAK below causes seg fault
  // free_match_collection(match_collection);

  carp(CARP_INFO, "crux-analyze-matches finished.");
  exit(0);
}

/*  ****************** Subroutines ****************/

/*
 * Outputs the matches in match_collection
 */
int output_matches(
    MATCH_COLLECTION_T* match_collection,
    SCORER_TYPE_T scorer_type
    ){
  // create match iterator, return match in sorted order of main_score type
  // TODO what is TRUE below?
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator(match_collection, scorer_type, TRUE);
  
  // print only up to max_rank_result of the matches
  int max_rank_result = get_int_parameter("max-sqt-result");

  // iterate over matches
  int match_count = 0;
  MATCH_T* match = NULL;
  while(match_iterator_has_next(match_iterator)){
    ++match_count;

    //// set max number of final scoring matches to print as output

    if(match_count >= max_rank_result){
      break;
    }
    
    match = match_iterator_next(match_iterator);
    // TODO what is TRUE below?
    print_match(match, stdout, TRUE, scorer_type);
  }
  return 0;
}

/**
 * Compare doubles
 */
int compare_doubles_descending(
    const void *a,
    const void *b
    ){
  double temp = *((double *)a) - *((double *)b);
  if (temp > 0){
    return -1;
  } else if (temp < 0){
    return 1;
  } else {
    return 0;
  }
}

/**
 * Perform Benjamini-Hochberg qvalue calculations on p-values generated
 * as in Klammer et al. (In Press) for PSMs in psm_result_folder, searched
 * against the sequence database in fasta_file. Requires that the match 
 * collection objects in the psm_result_folder have been scored using 
 * the p-value method (for now, only LOGP_BONF_WEIBULL_XCORR). 
 * There should be no decoy data sets in the directory.
 * \returns a MATCH_COLLECTION object
 */
MATCH_COLLECTION_T* run_nothing(
  char* psm_result_folder, 
  char* fasta_file,
  char* feature_file
  ){
  
  // create MATCH_COLLECTION_ITERATOR_T object
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_T* match = NULL;
  unsigned int number_features = 20;
  double* features = NULL;    

  FILE* feature_fh = NULL;
  
  // optional feature_file
  if (feature_file != NULL){
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
      return NULL;
    }
  }

  while(match_collection_iterator_has_next(match_collection_iterator)){

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    if (feature_fh == NULL){
      // Don't extract the features. Just return the match_collection
      return match_collection;
    }

    // create iterator, to register each PSM feature to Percolator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);

      // get the percolator features
      features = get_match_percolator_features(match, match_collection);
        
      fprintf(feature_fh, "%i\t",
          get_spectrum_first_scan(get_match_spectrum(match))
          );
      if (get_match_null_peptide(match) == FALSE){
        fprintf(feature_fh, "1\t%s\t", get_match_sequence(match));
      } else { 
        fprintf(feature_fh, "-1\t\t");
      };

      unsigned int feature_idx;
      for (feature_idx = 0; feature_idx < number_features; feature_idx++){
        if (feature_idx < number_features - 1){
          fprintf(feature_fh, "%.4f\t", features[feature_idx]);
        } else {
          fprintf(feature_fh, "%.4f\n", features[feature_idx]);
        }
      }
    }
    free_match_iterator(match_iterator);
  }
  fclose(feature_fh);

  free_match_collection_iterator(match_collection_iterator);
  return match_collection;
}

/**
 * Perform Benjamini-Hochberg qvalue calculations on p-values generated
 * as in Klammer et al. (In Press) for PSMs in psm_result_folder, searched
 * against the sequence database in fasta_file. Requires that the match 
 * collection objects in the psm_result_folder have been scored using 
 * the p-value method (for now, only LOGP_BONF_WEIBULL_XCORR). 
 * There should be no decoy data sets in the directory.
 * \returns a MATCH_COLLECTION object
 */
MATCH_COLLECTION_T* run_qvalue(
  char* psm_result_folder, 
  char* fasta_file 
  ){

  // double pi0 = get_double_parameter("pi0");
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_T* match = NULL;

  // array to store out pvalues
  const int length = MAX_PSMS;
  double* pvalues = (double*) malloc(sizeof(double) * length);
  int num_psms = 0;
  
  // create MATCH_COLLECTION_ITERATOR_T object
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);

  while(match_collection_iterator_has_next(match_collection_iterator)){

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    // create iterator, to register each PSM feature to Percolator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      pvalues[num_psms++] =  get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      if (num_psms >= MAX_PSMS){
        die("Too many psms in directory %s", psm_result_folder);
      }
    }

    // ok free & update for net set
    free_match_iterator(match_iterator);
  }

  free_match_collection_iterator(match_collection_iterator);

  // sort the - log p-values in descending order
  qsort(pvalues, num_psms, sizeof(double), compare_doubles_descending);

  double* qvalues = (double*) malloc(sizeof(double) * length);

  // work in negative log space, since that is where p- and qvalues end up
  double log_num_psms = - log(num_psms);
  double log_pi_0 = - log(get_double_parameter("pi0"));

  // convert the p-values into FDRs using Benjamini-Hochberg
  int idx;
  for (idx=0; idx < num_psms; idx++){
    carp(CARP_DETAILED_DEBUG, "pvalue[%i] = %.10f", idx, pvalues[idx]);
    int pvalue_idx = idx + 1; // start counting pvalues at 1
    double log_pvalue = pvalues[idx];

    double log_qvalue = 
      log_pvalue + log_num_psms - (-log(pvalue_idx)) + log_pi_0;
    qvalues[idx] = log_qvalue;
    carp(CARP_DETAILED_DEBUG, "no max qvalue[%i] = %.10f", idx, qvalues[idx]);
  }

  // convert the FDRs into q-values
  double max_log_qvalue = - BILLION;
  for (idx=num_psms-1; idx >= 0; idx--){
    if (qvalues[idx] > max_log_qvalue){
      max_log_qvalue = qvalues[idx];
    } else { // current q-value is <= max q-value
      // set to max q-value so far
      qvalues[idx] = max_log_qvalue; 
    }
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.10f", idx, qvalues[idx]);
  }

  // Iterate over the matches again
  match_collection_iterator = 
     new_match_collection_iterator(psm_result_folder, fasta_file);

  while(match_collection_iterator_has_next(match_collection_iterator)){

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    // create iterator, to register each PSM feature to Percolator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);

    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      double log_pvalue = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      carp(CARP_DETAILED_DEBUG, "- log pvalue  = %.6f", log_pvalue);
      
      // get the index of the p-value in the sorted list
      // FIXME slow, but it probably doesn't matter
      int idx;
      int pvalue_idx = -1;
      for (idx=0; idx < num_psms; idx++){
        double element = pvalues[idx];
        if ((element - EPSILON <= log_pvalue) &&
            (element + EPSILON >= log_pvalue)){
          pvalue_idx = idx; // start counting pvalues at 1
          break;
        }
      }
      
      set_match_score(match, LOGP_QVALUE_WEIBULL_XCORR, qvalues[pvalue_idx]);
    }

    // ok free & update for net set
    free_match_iterator(match_iterator);
    break; // just do the first match collection, which is the target matches
  }

   set_match_collection_scored_type(match_collection, 
       LOGP_QVALUE_WEIBULL_XCORR, TRUE);

  // free the match iterator, return match in sorted order of main_score type
  free_match_collection_iterator(match_collection_iterator);
  free(pvalues);
  free(qvalues);

  return match_collection;
}


/**
 * Runs Lukas Kall's percolator on the PSMs in the psm_result_folder
 * for a search against the sequence database fasta_file. Optionally 
 * puts the percolator PSM feature vectors into feature_file, if it is 
 * not NULL
 * \returns a MATCH_COLLECTION object
 */
MATCH_COLLECTION_T* run_percolator(
  char* psm_result_folder, 
  char* fasta_file, 
  char* feature_file){ 

  ALGORITHM_TYPE_T algorithm = PERCOLATOR_ALGORITHM;
  unsigned int number_features = 20;
  double* features = NULL;    
  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi0");
  char** feature_names = generate_feature_name_array(algorithm);
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_COLLECTION_T* target_match_collection = NULL;
  MATCH_T* match = NULL;
  FILE* feature_fh = NULL;
  int set_idx = 0;
  
  // optional feature_file
  if (feature_file != NULL){
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
      return NULL;
    }
  }

  carp(CARP_DETAILED_DEBUG, "Created feature file");

  // create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);

  carp(CARP_DETAILED_DEBUG, "Created the match collection iterator");

  // iterate over each, TARGET, DECOY 1..3 match_collection sets
  int iterations = 0;
  while(match_collection_iterator_has_next(match_collection_iterator)){
    
    carp(CARP_DEBUG, "Match collection iteration: %i" , iterations++);

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);
    
    // intialize percolator, using information from first match_collection
    if(set_idx == 0){
      // the first match_collection is the target_match_collection
      target_match_collection = match_collection;

      // result array that stores the algorithm scores
      results_q = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      results_score = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      
      // Call that initiates percolator
      pcInitiate(
          (NSet)get_match_collection_iterator_number_collections(
                  match_collection_iterator), 
          number_features, 
          get_match_collection_match_total(match_collection), 
          feature_names, 
          pi0);
      
      // Call that sets verbosity level
      // 0 is quiet, 2 is default, 5 is more than you want
      if(verbosity < CARP_ERROR){
        pcSetVerbosity(0);
      }    
      else if(verbosity < CARP_INFO){
        pcSetVerbosity(1);
      }
      else{
        pcSetVerbosity(5);
      }
    }

    // create iterator, to register each PSM feature to Percolator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      // Register PSM with features to Percolator    
      features = get_match_percolator_features(match, match_collection);

      if (feature_fh != NULL){
        
        fprintf(feature_fh, "%i\t",
            get_spectrum_first_scan(get_match_spectrum(match))
            );
        if (get_match_null_peptide(match) == FALSE){
          fprintf(feature_fh, "1\t");
        } else { 
          fprintf(feature_fh, "-1\t");
        };

        unsigned int feature_idx;
        for (feature_idx = 0; feature_idx < number_features; feature_idx++){
          if (feature_idx < number_features - 1){
            fprintf(feature_fh, "%.4f\t", features[feature_idx]);
          } else {
            fprintf(feature_fh, "%.4f\n", features[feature_idx]);
          }
        }
      }
      
      pcRegisterPSM((SetType)set_idx, 
                    NULL, // no sequence used
                    features);
      
      free(features);
    }

    // ok free & update for net set
    // MEMLEAK 
    free_match_iterator(match_iterator);

    // don't free the target_match_collection
    if(set_idx != 0){
      free_match_collection(match_collection);
    }

    ++set_idx;
  } // end iteratation over each, TARGET, DECOY 1..3 match_collection sets

  if (feature_fh != NULL){
    fclose(feature_fh);
  }
  
  /***** PERCOLATOR run *********/

  // Start processing
  pcExecute(); 
  
  /** Retrieving target scores and qvalues after 
   *  processing, the array should be numSpectra long and will be filled in 
   *  the same order as the features were inserted */
  pcGetScores(results_score, results_q); 
       
  // fill results for Q_VALUE
  fill_result_to_match_collection(
      target_match_collection, results_q, Q_VALUE, TRUE);
  
  // fill results for PERCOLATOR_SCORE
  fill_result_to_match_collection(
      target_match_collection, results_score, PERCOLATOR_SCORE, FALSE);
   
  // Function that should be called after processing finished
  pcCleanUp();
  
  // TODO put free back in. took out because claimed it was double free
  // free names
  /*unsigned int name_idx;
  for(name_idx=0; name_idx < number_features; ++name_idx){
    free(feature_names[name_idx]);
  }
  free(feature_names);
  */

  free(results_q);
  free(results_score);
  free_match_collection_iterator(match_collection_iterator);

  // TODO put free back in. took out because glibc claimed it was corrupted
  // double linked list
  // free_parameters();
  return target_match_collection;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
