/*****************************************************************************
 * \file match_analysis.c
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an optional parameter file, 
 * search all the spectrum against the peptides in the sequence database, and return high scoring peptides. 
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

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){
  char* usage = parse_arguments_get_usage("search_spectra");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  //print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  //FIXME uncomment this print if want to print usage whenever error message is printed
  //fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

int main(int argc, char** argv){
  //optional
  int verbosity = CARP_ERROR;
  char* parameter_file = NULL;
  char* psm_algorithm = "percolator";
  char* psm_result_folder = NULL;
  
  //required
  char* fasta_file = NULL;

  //parsing variables
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
    "The analysis algorithm to use. percolator|retention-czar|all",
    (void *) &psm_algorithm,
    STRING_ARG); 

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG); 
  
  parse_arguments_set_opt(
    "match-output-folder", 
    "The name of folder in which all the psm result files are located.",
    (void *) &psm_result_folder, 
    STRING_ARG);
  
  /* Define required command line arguments */
  parse_arguments_set_req(
    "fasta-file", 
    "The name of the file (in fasta format) from which to retrieve proteins and peptides.",
    (void *) &fasta_file, 
    STRING_ARG);
  
  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    ALGORITHM_TYPE_T algorithm = PERCOLATOR;
    MATCH_T* match = NULL;
    unsigned int name_idx = 0;
    
    //set max number of final scoring matches to print as output
    int max_rank_result = get_int_parameter("max-rank-result", 500);
    
    //set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity", "verbosity level must be between 0-100");
    }
    
    //set verbosity
    set_verbosity_level(verbosity);

    //parse and update parameters
    parse_update_parameters(parameter_file);
    
    //select algorithm
    if(strcmp(get_string_parameter_pointer("algorithm"), "percolator")== 0){
      algorithm = PERCOLATOR;
    }
    else if(strcmp(get_string_parameter_pointer("algorithm"), "retention-czar")== 0){
      algorithm = CZAR;
    }
    else if(strcmp(get_string_parameter_pointer("algorithm"), "all")== 0){
      algorithm = ALL;
    }
    else{
      wrong_command(psm_algorithm, "The analysis algorithm to use. percolator|retention-czar|all");
    }

    //parameters are now confirmed, can't be changed
    parameters_confirmed();
    
    unsigned int number_features = 20;
    double* results_q = NULL;
    double* results_score = NULL;
    double pi0 = 0.9; ///???
    NSet sets  = TWO_SETS;
    char** feature_names = generate_feature_name_array(algorithm);
    double* features = NULL;
    char* sequence = NULL;

    //create MATCH_COLLECTION_T object and read in the output PSM results.
    MATCH_COLLECTION_T* match_collection = new_match_collection_psm_output(psm_result_folder, fasta_file);
    //result array that stores the algorithm scores
    results_q = (double*)mycalloc(get_match_collection_match_total(match_collection), sizeof(double));
    results_score = (double*)mycalloc(get_match_collection_match_total(match_collection), sizeof(double));
    
    //Process PSM output and generate the run specific features
    //process_run_specific_features(match_collection);
    
    
    // Call that initiates percolator
    pcInitiate(sets, number_features, get_match_collection_match_total(match_collection), feature_names, pi0);
    

    //Call that sets verbosity level
    //0 is quiet, 2 is default, 5 is more than you want
    if(verbosity < CARP_ERROR){
      pcSetVerbosity(0);
    }    
    else if(verbosity < CARP_INFO){
      pcSetVerbosity(1);
    }
    else{
      pcSetVerbosity(2);
    }

    //create iterator, to register each PSM feature to Percolator
    MATCH_ITERATOR_T* match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      
      //Register PSM with features to Percolator    
      features = get_match_percolator_features(match, match_collection);
      sequence = get_match_sequence(match);
      
      pcRegisterPSM(TARGET, 
                    sequence, 
                    features);
     
      pcRegisterPSM(DECOY1, 
                    sequence, 
                    features);
      
      free(sequence);
      free(features);
    }

    /***** PERCOLATOR run *********/

    // Function called when we want to start processing
    pcExecute(); 
    
    /****************************/

    free_match_iterator(match_iterator);

    //create new iterator, to print result in sorted order of Q_Value
    //match_iterator = new_match_iterator(match_collection, Q_VALUE, TRUE);
    
    
    /** Function called when retrieving target scores and q-values after processing,
     * the array should be numSpectra long and will be filled in the same order
     * as the features were inserted */
    pcGetScores(results_score, results_q); 
    
    //fill results for Q_VALUE
    fill_result_to_match_collection(match_collection, results_q, Q_VALUE);
    
    //fill results for PERCOLATOR_SCORE
    fill_result_to_match_collection(match_collection, results_score, PERCOLATOR_SCORE);
    
    //create match iterator, TRUE: return match in sorted order of main_score type
    match_iterator = new_match_iterator(match_collection, Q_VALUE, TRUE);
    
    //iterate over matches
    int match_count = 0;
    while(match_iterator_has_next(match_iterator)){
      ++match_count;
      match = match_iterator_next(match_iterator);
      print_match(match, stdout, TRUE, Q_VALUE);
      
      //print only up to max_rank_result of the matches
      if(match_count >= max_rank_result){
        break;
      }
    }
    
    // Function that should be called after processing finished
    pcCleanUp();
    
    //free names
    for(; name_idx < number_features; ++name_idx){
      free(feature_names[name_idx]);
    }
    free(feature_names);

    free(results_q);
    free(results_score);
    free_match_iterator(match_iterator);
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
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
