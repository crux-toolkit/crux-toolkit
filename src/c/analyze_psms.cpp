/**
 * \file analyze_psms.cpp
 */

#include "analyze_psms.h"
#include "ComputeQValues.h"
#include "QRanker.h"
#include "Percolator.h"
#include "MatchCollectionIterator.h"
#include "MatchIterator.h"


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

  CruxApplication* application = NULL;

  // Do some common initialization stuff.
  switch(command) {
  case QVALUE_COMMAND:
  {
    application = new ComputeQValues();
    application->initialize(argument_list, num_arguments,
      qvalue_option_list, qvalue_num_options, argc, argv);
  }
    break;
  case PERCOLATOR_COMMAND:
  {
    application = new Percolator();
    application->initialize(argument_list, num_arguments,
      percolator_option_list, percolator_num_options, argc, argv);
  }
    break;
  case QRANKER_COMMAND:
  {
    application = new QRanker();
    application->initialize(argument_list, num_arguments,
                   qranker_option_list, qranker_num_options, argc, argv);
  }
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  // Get required arguments
  char* protein_database_name = get_string_parameter("protein database");
  char* input_directory = get_string_parameter("search results directory");

  // Prepare the output files.
  OutputFiles output(application);

  // Perform the analysis.
  MatchCollection* match_collection = NULL;
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

  delete match_collection;
  output.writeFooters();
  // clean up
  free(input_directory);
  free(protein_database_name);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  string name = application->getName();

  carp(CARP_INFO, "Finished crux %s.", name.c_str());
  delete application;
}

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
MatchCollection* run_percolator_or_qranker(
  COMMAND_T command,                                          
  char* input_directory, 
  char* fasta_file, // actually fasta or index
  OutputFiles& output){ 

  double* features = NULL;    
  double* results_q = NULL;
  double* results_score = NULL;
  double pi0 = get_double_parameter("pi-zero");
  char** feature_names = generate_feature_name_array();
  MatchIterator* match_iterator = NULL;
  MatchCollection* match_collection = NULL;
  MatchCollection* target_match_collection = NULL;
  Match* match = NULL;
  int set_idx = 0;
  
  output.writeFeatureHeader(feature_names, NUM_FEATURES);

  // Create MATCH_COLLECTION_ITERATOR_T object
  // which will read in the serialized output PSM results and return
  // first the match_collection of TARGET followed by 
  // the DECOY* match_collections.
  int num_decoys = 0;
  MatchCollectionIterator* match_collection_iterator =
    new MatchCollectionIterator(input_directory, fasta_file, &num_decoys);
  // Create an array with counts of spectra in each match collection.
  int num_sets = match_collection_iterator->getNumberCollections();
  int* num_spectra = (int*)mycalloc(num_sets, sizeof(int));
  int iterations = 0;
  while(match_collection_iterator->hasNext()){
    match_collection = 
      match_collection_iterator->next();
    num_spectra[iterations] = 
      match_collection->getMatchTotal();
    iterations++;
  }

  // get from the match files the columns to print in the output files
  const vector<bool>& cols_to_print = 
    match_collection_iterator->getColsInFile();
  output.writeHeaders(cols_to_print);

  // Reset the iterator. (FIXME: There should be a function to do this!)
  delete match_collection_iterator;
  num_decoys = 0;
  match_collection_iterator =
    new MatchCollectionIterator(input_directory, fasta_file, &num_decoys);

  // iterate over each, TARGET, DECOY 1..3 match_collection sets
  iterations = 0;
  while(match_collection_iterator->hasNext()){
    
    carp(CARP_DEBUG, "Match collection iteration: %i" , iterations++);

    // get the next match_collection
    match_collection = match_collection_iterator->next();
    
    // intialize percolator, using information from first match_collection
    if (set_idx == 0) {
      // the first match_collection is the target_match_collection
      target_match_collection = match_collection;

      // result array that stores the algorithm scores
      results_q = (double*)mycalloc(
          match_collection->getMatchTotal(), sizeof(double));
      results_score = (double*)mycalloc(
          match_collection->getMatchTotal(), sizeof(double));
      
      // Call that initiates q-ranker or percolator.
      switch (command) {
      case PERCOLATOR_COMMAND:
        pcInitiate(
          (NSet)match_collection_iterator->getNumberCollections(), 
          NUM_FEATURES, 
          num_spectra, 
          feature_names, 
          pi0);
        break;
      case QRANKER_COMMAND:
        qcInitiate(
            (NSet)match_collection_iterator->getNumberCollections(),
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
      case MISC_COMMAND:
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
    carp(CARP_DETAILED_DEBUG, "Registering PSMs");
    match_iterator = new MatchIterator(match_collection, XCORR, false);
    
    while(match_iterator->hasNext()){
      match = match_iterator->next();
      // Register PSM with features
      features = match->getPercolatorFeatures(match_collection);

      output.writeMatchFeatures(match, features, NUM_FEATURES);
      switch (command) {
      case PERCOLATOR_COMMAND:
        pcRegisterPSM((SetType)set_idx, 
                      NULL, // no sequence used
                      features);
        break;
      case QRANKER_COMMAND:
        qcRegisterPSM((SetType)set_idx,
                      match->getSequenceSqt(),
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
    delete match_iterator;

    // don't free the target_match_collection
    if(set_idx != 0){
      delete match_collection;
    }

    ++set_idx;
  } // end iteratation over each, TARGET, DECOY 1..3 match_collection sets

  carp(CARP_DETAILED_DEBUG, "Processing PSMs");
  // Start processing
  switch (command) {
  case PERCOLATOR_COMMAND:
    pcExecute(); 
    pcGetScores(results_score, results_q); 
    target_match_collection->fillResult(
      results_q, PERCOLATOR_QVALUE, true);
    target_match_collection->fillResult(
      results_score, PERCOLATOR_SCORE, false);
    pcCleanUp();
    break;
  case QRANKER_COMMAND:
    qcExecute(!get_boolean_parameter("no-xval")); 
    qcGetScores(results_score, results_q); 
    target_match_collection->fillResult(
        results_q, QRANKER_QVALUE, true);
    target_match_collection->fillResult(
        results_score, QRANKER_SCORE, false);
    qcCleanUp();
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  carp(CARP_DETAILED_DEBUG, "Writing matches");
  output.writeMatches(target_match_collection);

  carp(CARP_DETAILED_DEBUG, "analyze_psms: cleanup");

  // free names
  unsigned int name_idx;
  for(name_idx=0; name_idx < NUM_FEATURES; ++name_idx){
    free(feature_names[name_idx]);
  }
  free(feature_names);
  free(results_q);
  free(results_score);
  delete match_collection_iterator;

  // TODO put free back in. took out because glibc claimed it was corrupted
  // double linked list
  // free_parameters();

  return target_match_collection;
}





