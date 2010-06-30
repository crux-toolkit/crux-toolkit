/**
 * \file crux-main.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: November 24, 2008
 * \brief The starting point for the main crux program.
 *
 * Usage is "crux [command] [options] [arguments]" where command
 * is one of the primary crux commands.
 **/

#include "crux-main.h"
#include "crux-utils.h" // Need to get definition of NUM_FEATURES.

const char* usage_str = "Usage: crux <command> [options] <argument>\n"
"\n"
"Crux supports the following commands:\n"
"  create-index        Create an index for all peptides in a fasta file.\n"
"  search-for-matches  Search a collection of spectra against a sequence\n"
"                      database, returning a collection of peptide-spectrum\n"
"                      matches (PSMs) scored by XCorr.\n"
"  sequest-search      Similar to search-for-matches but use Sp as a \n"
"                      preliminary score followed by XCorr.\n"
"  compute-q-values    Assign a q-value, which is a statistical confidence\n"
"                      measure that accounts for multiple testing, to each\n"
"                      PSM in a given set.\n" 
"  percolator          Analyze a collection of PSMs to target and decoy\n"
"                      sequences using the percolator algorithm.\n"
"  q-ranker            Analyze a collection of PSMs using the Q-ranker\n"
"                      algorithm.\n"
"  print-processed-spectra\n"
"                      Write a new ms2 file with all of the same spectra\n"
"                      with only the peaks used for computing xcorr.\n"
"  search-for-xlinks   Search a collection of spectra against a sequence\n"
"                      database, returning a collection of matches\n"
"                      corresponding to linear and cross-linked peptides\n"
"                      scored by XCorr.\n"
"  version             Print the Crux version number to standard output,\n"
"                      then exit.\n"
"\n"
"Options and arguments are specific to each command. Type 'crux <command>'\n"
"for details.\n"
; 

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
static MATCH_COLLECTION_T* run_percolator_or_qranker(
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
  // N.B. This array is actually only needed by q-ranker.
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
          get_match_collection_match_total(match_collection), 
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
      default:
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
        target_match_collection, results_q, Q_VALUE, TRUE);
    fill_result_to_match_collection(
        target_match_collection, results_score, PERCOLATOR_SCORE, FALSE);
    pcCleanUp();
    break;
  case QRANKER_COMMAND:
    qcExecute(!get_boolean_parameter("no-xval")); 
    qcGetScores(results_score, results_q); 
    fill_result_to_match_collection(
        target_match_collection, results_q, QRANKER_Q_VALUE, TRUE);
    fill_result_to_match_collection(
        target_match_collection, results_score, QRANKER_SCORE, FALSE);
    qcCleanUp();
    break;
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

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


/**
 * \brief Takes a directory containing PSM files and a protein index
 * and analyzes the PSMs using compute-q-values, percolator or q-ranker.
 */
static void analyze_matches_main(
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
  output.writeHeaders();

  // Perform the analysis.
  MATCH_COLLECTION_T* match_collection = NULL;
  switch(command) {
  case QVALUE_COMMAND:
    match_collection = run_qvalue(input_directory,
				  protein_database_name);
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

  carp(CARP_INFO, "Outputting matches.");
  output.writeMatches(match_collection);

  // MEMLEAK below causes seg fault (or used to)
  // free_match_collection(match_collection);

  // clean up
  free(input_directory);
  free(protein_database_name);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  char* name = command_type_to_command_line_string(command);
  carp(CARP_INFO, "Finished crux %s.", name);
}


int main(int argc, char** argv){

  // check the syntax for crux <operation>
  if( argc < 2 ){
    carp(CARP_FATAL, usage_str);
  }

  // determine the operation
  char* op_string = argv[1];
  COMMAND_T command = string_to_command_type(op_string);

  // call the appropriate function 
  // passing the command line minus the first token ('crux')
  switch(command){
  case INDEX_COMMAND:
    create_index_main(argc-1, argv+1);
    break;

  case SEARCH_COMMAND:
    search_main(argc-1, argv+1);
    break;

  case SEQUEST_COMMAND:
    sequest_search_main(argc-1, argv+1);
    break;

  case PROCESS_SPEC_COMMAND:
    print_processed_spectra_main(argc-1, argv+1);
    break;
  
  case XLINK_SEARCH_COMMAND:
    xlink_search_main(argc-1, argv+1);
    break;

  case QVALUE_COMMAND:
  case QRANKER_COMMAND:
  case PERCOLATOR_COMMAND:
    analyze_matches_main(command, argc-1, argv+1);
    break;

  case VERSION_COMMAND:
    printf("Crux version %s\n", VERSION);
    break;    

  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command '%s'\n%s", op_string, usage_str);
    break;

  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;

  }

  exit (0);
}// end main















