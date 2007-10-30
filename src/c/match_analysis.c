/*****************************************************************************
 * \file match_analysis.c
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an 
 *              optional parameter file, search all the spectrum against 
 *              the peptides in the sequence database, and return the high 
 *              scoring peptides. 
 * REVISION: 
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"
#include "peptide.h"
#include "protein.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "crux-utils.h"
#include "scorer.h"
#include "objects.h"
#include "parameter.h"
#include "match.h"
#include "match_collection.h"
#include "PercolatorCInterface.h"

#define MAX_PSMS 10000000
#define EPSILON 0.000001

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){
  char* usage = parse_arguments_get_usage("search_spectra");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  // print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

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
  int max_rank_result = get_int_parameter("max-rank-result");

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

/***********************************************************************/
int main(int argc, char** argv){
  // optional
  int verbosity = CARP_ERROR;
  char* parameter_file = "crux.params";
  char* psm_algorithm = "percolator";
  char* psm_result_folder = NULL;
  
  // required
  char* fasta_file   = NULL;
  char* feature_file = NULL;

  // parsing variables
  int result = 0;
  char* error_message;
  
  /* Define optional command line arguments */   
  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);
  
  parse_arguments_set_opt(
    "algorithm",
    "The analysis algorithm to use. percolator|retention-czar|qvalue|all",
    (void *) &psm_algorithm,
    STRING_ARG); 

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG); 
  
  parse_arguments_set_req(
    "match-output-folder", 
    "The name of folder in which all the psm result files are located.",
    (void *) &psm_result_folder, 
    STRING_ARG);
  
  /* Define required command line arguments */
  parse_arguments_set_req(
    "fasta-file", 
    "The name of the file (in fasta format) from which to retrieve proteins "
    "and peptides.",
    (void *) &fasta_file, 
    STRING_ARG);
  
  parse_arguments_set_opt(
    "feature-file",
    "Optional file in which to write the features",
    (void *) &feature_file,
    STRING_ARG); 

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
        
    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity", "verbosity level must be between 0-100");
    }
    
    // set verbosity
    set_verbosity_level(verbosity);

    // parse and update parameters
    parse_update_parameters(parameter_file);
    
    // always use index for match_analysis!
    set_string_parameter("use-index", "T");
    
    // parameters are now confirmed, can't be changed
    parameters_confirmed();
    
    ALGORITHM_TYPE_T algorithm = PERCOLATOR;

    // select algorithm TODO put this in a routine in parameter.c
    char* algorithm_string = get_string_parameter_pointer("algorithm");
    if(strcmp(algorithm_string, "percolator")== 0){
      algorithm = PERCOLATOR;
    }
    else if(strcmp(algorithm_string, "retention-czar")== 0){
      algorithm = CZAR;
    }
    else if(strcmp(algorithm_string, "qvalue")== 0){
      algorithm = QVALUE;
    }
    else if(strcmp(algorithm_string, "all")== 0){
      algorithm = ALL;
    }
    else{
      wrong_command(psm_algorithm, 
        "The analysis algorithm to use. percolator|retention-czar|all");
    }

    MATCH_COLLECTION_T* match_collection = NULL;
    SCORER_TYPE_T scorer_type = 0;
    if (algorithm == PERCOLATOR){
      carp(CARP_INFO, "Running percolator");
      match_collection = run_percolator(
          psm_result_folder, fasta_file, feature_file); 
      scorer_type = PERCOLATOR_SCORE;
    } else if (algorithm == QVALUE){
      carp(CARP_INFO, "Running q-value");
      match_collection = run_qvalue(psm_result_folder, fasta_file);
      scorer_type = LOGP_QVALUE_WEIBULL_XCORR;
    }
    output_matches(match_collection, scorer_type);
    free_match_collection(match_collection);
   }
  else{
    char* usage = parse_arguments_get_usage("match_analysis");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
  }
  exit(0);
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
    carp(CARP_DETAILED_DEBUG, "pvalue[%i] = %.6f", idx, pvalues[idx]);
    int pvalue_idx = idx + 1; // start counting pvalues at 1
    double log_pvalue = pvalues[idx];

    double log_qvalue = 
      log_pvalue + log_num_psms - (-log(pvalue_idx)) + log_pi_0;
    qvalues[idx] = log_qvalue;
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
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.6f", idx, qvalues[idx]);
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
      int pvalue_idx;
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
  // MEMLEAK this is free'ing the database, however, we still have pointers
  // to the db in the protein objects. Not sure this how we want to deal
  // with the database, however, since it creates circular references
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

  ALGORITHM_TYPE_T algorithm = PERCOLATOR;
  unsigned int number_features = 20;
  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi0");
  char** feature_names = generate_feature_name_array(algorithm);
  double* features = NULL;    
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

  // create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);

  // iterate over each, TARGET, DECOY 1..3 match_collection sets
  int iterations = 0;
  while(match_collection_iterator_has_next(match_collection_iterator)){
    
    carp(CARP_INFO, "Match collection iteration: %i" , iterations++);

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
  free_match_iterator(match_iterator);

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
