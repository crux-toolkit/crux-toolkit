/*************************************************************************//**
 * \file match_analysis.c
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * \brief  Given as input a directory containing binary psm files,
 * a protein database, and an optional parameter file analyze the
 * matches (with percolator or q-value) and return scores indicating
 * how good the matches are. 
 *
 * Handles at most 4 files (target and decoy).  Expects psm files to
 * end with the extension '.csm' and decoys to end with
 * '-decoy#.csm'.  Multiple target files in the given directory are
 * concatinated together and presumed to be non-overlaping parts of
 * the same ms2 file. 
 * 
 * $Revision: 1.51 $
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
#include "parse_arguments.h" 
#include "spectrum_collection.h"
#include "generate_peptides_iterator.h"
#include "scorer.h"
#include "match.h"
#include "match_collection.h"
#include "hit.h"
#include "hit_collection.h"
#include "PercolatorCInterface.h"

#define MAX_PSMS 10000000
// 14th decimal place
#define EPSILON 0.00000000000001 
#define NUM_ANALYSIS_OPTIONS 7
#define NUM_ANALYSIS_ARGUMENTS 2

/* 
 * Private function declarations.  Details below
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

void print_sqt_file(
  MATCH_COLLECTION_T* match_collection,
  SCORER_TYPE_T scorer_type,
  SCORER_TYPE_T second_scorer_type
  );

  
/**
 * \brief crux-analyze-matches: takes in a directory containing binary
 * psm files and a protein index and analyzes the psms.
 */
int main(int argc, char** argv){

  /* Define command line arguments */
  int num_options = NUM_ANALYSIS_OPTIONS;
  char* option_list[NUM_ANALYSIS_OPTIONS] = {
    "version",
    "verbosity",
    "parameter-file",
    "algorithm",
    "feature-file",
    "overwrite"
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
  //set_verbosity_level(get_int_parameter("verbosity"));

  /* Get arguments */
  //  char* psm_file = get_string_parameter("psm file");
  char* psm_file = get_string_parameter("psm-folder");
  char* fasta_file = get_string_parameter("protein input");//rename
  char* feature_file = get_string_parameter("feature-file");

  /* Get options */
  ALGORITHM_TYPE_T algorithm_type = get_algorithm_type_parameter("algorithm");
  SCORER_TYPE_T scorer_type = PERCOLATOR_SCORE;
  SCORER_TYPE_T second_scorer_type = Q_VALUE;
  MATCH_COLLECTION_T* match_collection = NULL;

  /* Perform the analysis */
  switch(algorithm_type){
  case PERCOLATOR_ALGORITHM:
    carp(CARP_INFO, "Running percolator");
    match_collection = run_percolator(psm_file,
                                      fasta_file,
                                      feature_file);
    scorer_type = PERCOLATOR_SCORE;
    second_scorer_type = Q_VALUE;
    break;
    
  case QVALUE_ALGORITHM:
    carp(CARP_INFO, "Running qvalue");
    match_collection = run_qvalue(psm_file, fasta_file);
    scorer_type =  LOGP_QVALUE_WEIBULL_XCORR; 
    second_scorer_type = XCORR; // could it be other?
    break;
    
  case NO_ALGORITHM:
    carp(CARP_INFO, "No analysis algorithm chosen.");
    match_collection = run_nothing(psm_file,
                                   fasta_file,
                                   feature_file);
    scorer_type = XCORR; // TODO put in something to default to the primary
    second_scorer_type = SP;
    // score in the run
    break;
  default:
    ;
  }  
  
  carp(CARP_INFO, "Outputting matches.");
  //  output_matches(match_collection, scorer_type);
  print_sqt_file(match_collection, scorer_type, second_scorer_type);

  // MEMLEAK below causes seg fault (or used to)
  // free_match_collection(match_collection);

  /*
   The method new_hit_collection_from_match_collection below is an
   example of how one might assemble peptide identifications (matches)
   into protein identifications (hits). 
   Unfortunately it doesn't work that well, so it's commented out. But
   some of the functionality one would need is hopefully there. 
   */
  /* carp(CARP_INFO, "Assembling matches into protein hits");
  HIT_COLLECTION_T* hit_collection 
    = new_hit_collection_from_match_collection(match_collection);
  carp(CARP_INFO, "Outputting protein hits");
  print_hit_collection(stdout, hit_collection);
  free_hit_collection(hit_collection);
  */

  // clean up
  free(psm_file);
  free(fasta_file);//rename
  free(feature_file);


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

    // set max number of final scoring matches to print as output

    if(match_count >= max_rank_result){
      break;
    }
    
    match = match_iterator_next(match_iterator);
    print_match(match, stdout, TRUE, // TRUE==print sequence
                scorer_type);
  }
  free_match_iterator(match_iterator);
  return 0;
}

/*
 */
void print_sqt_file(
  MATCH_COLLECTION_T* match_collection,
  SCORER_TYPE_T scorer,
  SCORER_TYPE_T second_scorer
  ){

  // get filename and open file
  char* sqt_filename = get_string_parameter("sqt-output-file");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  FILE* sqt_file = create_file_in_path( sqt_filename, NULL, overwrite );

  // print header
  int num_proteins = get_match_collection_num_proteins(match_collection);
  print_sqt_header( sqt_file, "target", num_proteins, TRUE);

  ALGORITHM_TYPE_T algorithm_type = get_algorithm_type_parameter("algorithm");
  char algorithm_str[64];
  algorithm_type_to_string(algorithm_type, algorithm_str);

  fprintf(sqt_file, "H\tComment\tmatches analyzed by %s\n", algorithm_str);

  // get match iterator sorted by spectrum
  MATCH_ITERATOR_T* match_iterator = 
    new_match_iterator_spectrum_sorted(match_collection, scorer);

  // print each spectrum only once, keep track of which was last printed
  int cur_spectrum_num = -1;
  int cur_charge = 0;
  int match_counter = 0;
  int max_matches = get_int_parameter("max-sqt-result");

  // for all matches
  while( match_iterator_has_next(match_iterator) ){

    // get match and spectrum
    MATCH_T* match = match_iterator_next(match_iterator);
    SPECTRUM_T* spectrum = get_match_spectrum(match);
    int this_spectrum_num = get_spectrum_first_scan(spectrum);
    int charge = get_match_charge(match);

    carp(CARP_DETAILED_DEBUG, 
         "SQT printing scan %i (current %i), charge %i (current %i)", 
         this_spectrum_num, cur_spectrum_num, charge, cur_charge);

    // if this spectrum has not been printed...
    if( cur_spectrum_num != this_spectrum_num
        || cur_charge != charge){

      carp(CARP_DETAILED_DEBUG, "Printing new S line");
      // print S line to sqt file
      cur_spectrum_num = this_spectrum_num;
      cur_charge = charge;
      int num_peptides = get_match_ln_experiment_size(match);
      num_peptides = expf(num_peptides);

      print_spectrum_sqt(spectrum, sqt_file, num_peptides, charge);

      // print match to sqt file
      print_match_sqt(match, sqt_file, scorer, second_scorer);
      match_counter = 1;
    }
    // if this spectrum has been printed
    else{  
      if( match_counter < max_matches ){
        print_match_sqt(match, sqt_file, scorer, second_scorer);
        match_counter++;
      }
    }

  }// next match
  free_match_iterator(match_iterator);
  free(sqt_filename);

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
  MATCH_COLLECTION_T* match_collection = NULL; // to return
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
    }
  }

  while(match_collection_iterator_has_next(match_collection_iterator)){

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);

    if (feature_fh == NULL){
      // Don't extract the features. Just return the (first) match_collection
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
    }// next match for this collection
    free_match_iterator(match_iterator);

  }// next match collection (target, decoy, decoy...)
  fclose(feature_fh);

  free_match_collection_iterator(match_collection_iterator);
  return match_collection;  // this is the last returned by the
                            // iterator, likely a decoy, yes?
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

    // create iterator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);

    // for each match, get p-value    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);

      FLOAT_T score = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      carp(CARP_DETAILED_DEBUG, "p-value is %f", score);
      if( score == P_VALUE_NA ){// ignore unscored psms
        continue;
      }
      //pvalues[num_psms++] =  get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      pvalues[num_psms++] =  score;
      if (num_psms >= MAX_PSMS){
        carp(CARP_ERROR, "Too many psms in directory %s", psm_result_folder);
      }
    }

    // ok free & update for next set
    free_match_iterator(match_iterator);
  }// next match collection
  carp(CARP_DEBUG, "Read in %i p-values", num_psms);

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

    // create iterator
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);

    // for each match, convert p-value to q-value
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      double log_pvalue = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
      carp(CARP_DETAILED_DEBUG, "- log pvalue  = %.6f", log_pvalue);

      // if p-value wasn't calculated, set q-value as nan
      if( log_pvalue == P_VALUE_NA ){
        set_match_score(match, LOGP_QVALUE_WEIBULL_XCORR, sqrt(-1) );
        continue;
      }
      
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

    // ok free & update for next set
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
 * \brief Analyze matches using the percolator algorithm
 * 
 * Runs Lukas Kall's percolator on the PSMs in the psm_result_folder
 * for a search against the sequence database fasta_file. Optionally 
 * puts the percolator PSM feature vectors into feature_file, if it is 
 * not NULL.
 * \returns a pointer to a MATCH_COLLECTION_T object
 * \callgraph
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
  if(feature_file != NULL){  
    if((feature_fh = fopen(feature_file, "w")) == NULL){
      carp(CARP_FATAL, "Problem opening output file %s", feature_file);
    }
  }

  carp(CARP_DETAILED_DEBUG, "Created feature file");

  // create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(psm_result_folder, fasta_file);

  if( match_collection_iterator == NULL ){
    carp(CARP_FATAL, "Failed to create a match collection iterator");
  }
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

    // ok free & update for next set
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
  
  /* Retrieving target scores and qvalues after 
   * processing, the array should be numSpectra long and will be filled in 
   * the same order as the features were inserted */
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
  unsigned int name_idx;
  for(name_idx=0; name_idx < number_features; ++name_idx){
    free(feature_names[name_idx]);
  }
  free(feature_names);
  

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
