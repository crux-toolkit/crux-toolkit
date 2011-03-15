/**
 * \file analyze_psms.cpp
 */

#include "analyze_psms.h"

/**
 * \brief Takes a directory containing PSM files and a protein index
 * and analyzes the PSMs using compute-q-values, percolator or q-ranker.
 */
void analyze_matches_main(
  COMMAND_T command,
  int argc,
  char** argv
){

  // Define optional command line arguments.
  const char* qvalue_option_list[] = {
    "pi-zero",
    "verbosity",
    "parameter-file",
    "overwrite",
    "output-dir",
    "fileroot"
  };
  int qvalue_num_options = sizeof(qvalue_option_list)/sizeof(char*);
  const char* percolator_option_list[] = {
    "pi-zero",
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite"
  };
  int percolator_num_options = sizeof(percolator_option_list)/sizeof(char*);
  const char* qranker_option_list[] = {
    "pi-zero",
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite",
  };
  int qranker_num_options = sizeof(qranker_option_list)/sizeof(char*);

  // Define required command line arguments.
  const char* argument_list[] = {
    "protein database",
    "search results directory"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  // Do some common initialization stuff.
  switch(command) {
  case QVALUE_COMMAND:
    initialize_run(command, argument_list, num_arguments,
                   qvalue_option_list, qvalue_num_options, argc, argv);
    break;
  case PERCOLATOR_COMMAND:
    initialize_run(command, argument_list, num_arguments,
                   percolator_option_list, percolator_num_options, argc, argv);
    break;
  case QRANKER_COMMAND:
    initialize_run(command, argument_list, num_arguments,
                   qranker_option_list, qranker_num_options, argc, argv);
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  // Get required arguments
  char* protein_database_name = get_string_parameter("protein database");
  char* input_directory = get_string_parameter("search results directory");

  // Prepare the output files.
  OutputFiles output(command);

  // Perform the analysis.
  MATCH_COLLECTION_T* match_collection = NULL;
  switch(command) {
  case QVALUE_COMMAND:
    match_collection = run_qvalue(input_directory,
                                  protein_database_name,
                                  output);
    break;
  case PERCOLATOR_COMMAND:
  case QRANKER_COMMAND:
    match_collection = run_percolator_or_qranker(command,
                                                 input_directory,
                                                 protein_database_name,
                                                 output);
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  // NOTE: return this step once the db/index is created and freed here
  //carp(CARP_INFO, "Outputting matches.");
  //output.writeMatches(match_collection);

  free_match_collection(match_collection);
  output.writeFooters();
  // clean up
  free(input_directory);
  free(protein_database_name);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  char* name = command_type_to_command_line_string(command);
  carp(CARP_INFO, "Finished crux %s.", name);
  free(name);
}

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
  OutputFiles& output){ 

  double* features = NULL;    
  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi-zero");
  char** feature_names = generate_feature_name_array();
  MATCH_ITERATOR_T* match_iterator = NULL;
  MATCH_COLLECTION_T* match_collection = NULL;
  MATCH_COLLECTION_T* target_match_collection = NULL;
  MATCH_T* match = NULL;
  int set_idx = 0;
  
  output.writeFeatureHeader(feature_names, NUM_FEATURES);

  // Create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  int num_decoys = 0;
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator =
    new_match_collection_iterator(input_directory, fasta_file, &num_decoys);

  // Create an array with counts of spectra in each match collection.
  int num_sets = get_match_collection_iterator_number_collections(
      match_collection_iterator);
  int* num_spectra = (int*)mycalloc(num_sets, sizeof(int));
  int iterations = 0;
  while(match_collection_iterator_has_next(match_collection_iterator)){
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);
    num_spectra[iterations] = 
      get_match_collection_match_total(match_collection);
    iterations++;
  }

  // get from the match files the columns to print in the output files
  const vector<bool>& cols_to_print = 
    get_match_collection_iterator_cols_in_file(match_collection_iterator);
  output.writeHeaders(cols_to_print);

  // Reset the iterator. (FIXME: There should be a function to do this!)
  free_match_collection_iterator(match_collection_iterator);
  num_decoys = 0;
  match_collection_iterator =
    new_match_collection_iterator(input_directory, fasta_file, &num_decoys);

  // iterate over each, TARGET, DECOY 1..3 match_collection sets
  iterations = 0;
  while(match_collection_iterator_has_next(match_collection_iterator)){
    
    carp(CARP_DEBUG, "Match collection iteration: %i" , iterations++);

    // get the next match_collection
    match_collection = 
      match_collection_iterator_next(match_collection_iterator);
    
    // intialize percolator, using information from first match_collection
    if (set_idx == 0) {
      // the first match_collection is the target_match_collection
      target_match_collection = match_collection;

      // result array that stores the algorithm scores
      results_q = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      results_score = (double*)mycalloc(
          get_match_collection_match_total(match_collection), sizeof(double));
      
      // Call that initiates q-ranker or percolator.
      switch (command) {
      case PERCOLATOR_COMMAND:
        pcInitiate(
          (NSet)get_match_collection_iterator_number_collections(
                  match_collection_iterator), 
          NUM_FEATURES, 
          num_spectra, 
          feature_names, 
          pi0);
        break;
      case QRANKER_COMMAND:
        qcInitiate(
            (NSet)get_match_collection_iterator_number_collections(
                     match_collection_iterator), 
            NUM_FEATURES, 
            num_spectra, 
            feature_names, 
            pi0);
        break;
      case INDEX_COMMAND:
      case SEARCH_COMMAND:
      case SEQUEST_COMMAND:
      case QVALUE_COMMAND:
      case SPECTRAL_COUNTS_COMMAND:
      case PROCESS_SPEC_COMMAND:
      case XLINK_SEARCH_COMMAND:
      case VERSION_COMMAND:
      case INVALID_COMMAND:
      case NUMBER_COMMAND_TYPES:
        carp(CARP_FATAL, "Unknown command type.");
        break;
      }
      free(num_spectra);

      // Call that sets verbosity level
      // 0 is quiet, 2 is default, 5 is more than you want
      switch (command) {
      case PERCOLATOR_COMMAND:
        if(verbosity < CARP_ERROR){
          pcSetVerbosity(0);
        }    
        else if(verbosity < CARP_INFO){
          pcSetVerbosity(1); // FIXME
        }
        else{
          pcSetVerbosity(5); // FIXME
        }
        break;
      case QRANKER_COMMAND:
        if(verbosity < CARP_ERROR){
          qcSetVerbosity(0);
        }    
        else if(verbosity < CARP_INFO){
          qcSetVerbosity(1);
        }
        else{
          qcSetVerbosity(5);
        }
        break;
      default:
        carp(CARP_FATAL, "Unknown command type.");
        break;
      }
    }

    // create iterator, to register each PSM feature.
    match_iterator = new_match_iterator(match_collection, XCORR, FALSE);
    
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      // Register PSM with features
      features = get_match_percolator_features(match, match_collection);

      output.writeMatchFeatures(match, features, NUM_FEATURES);
      switch (command) {
      case PERCOLATOR_COMMAND:
        pcRegisterPSM((SetType)set_idx, 
                      NULL, // no sequence used
                      features);
        break;
      case QRANKER_COMMAND:
        qcRegisterPSM((SetType)set_idx,
                      get_match_sequence_sqt(match),
                      features);
        break;
      default:
        carp(CARP_FATAL, "Unknown command type.");
        break;
      }
      
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

  // Start processing
  switch (command) {
  case PERCOLATOR_COMMAND:
    pcExecute(); 
    pcGetScores(results_score, results_q); 
    fill_result_to_match_collection(
        target_match_collection, results_q, PERCOLATOR_QVALUE, TRUE);
    fill_result_to_match_collection(
        target_match_collection, results_score, PERCOLATOR_SCORE, FALSE);
    pcCleanUp();
    break;
  case QRANKER_COMMAND:
    qcExecute(!get_boolean_parameter("no-xval")); 
    qcGetScores(results_score, results_q); 
    fill_result_to_match_collection(
        target_match_collection, results_q, QRANKER_QVALUE, TRUE);
    fill_result_to_match_collection(
        target_match_collection, results_score, QRANKER_SCORE, FALSE);
    qcCleanUp();
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  output.writeMatches(target_match_collection);

  // free names
  unsigned int name_idx;
  for(name_idx=0; name_idx < NUM_FEATURES; ++name_idx){
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





