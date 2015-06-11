/**
 * \file analyze_psms.cpp
 */

#include <vector>
#include <functional>
#ifdef _MSC_VER
#include <iterator>
#endif
#include "analyze_psms.h"
#include "ComputeQValues.h"
#include "QRanker.h"
//#include "Percolator.h"
#include "MatchCollectionIterator.h"
#include "MatchCollectionParser.h"
#include "MatchIterator.h"
//#include "PercolatorCInterface.h"
#include "PosteriorEstimator.h"

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
/*  const char* percolator_option_list[] = {
    "pi-zero",
    "verbosity",
    "parameter-file",
    "fileroot",
    "feature-file",
    "output-dir",
    "overwrite"
  };
  int percolator_num_options = sizeof(percolator_option_list)/sizeof(char*);
*/
  // Define required command line arguments.
  const char* argument_list[] = {
    "target input"
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
/*
  case PERCOLATOR_COMMAND:
  {
    application = new Percolator();
    application->initialize(argument_list, num_arguments,
      percolator_option_list, percolator_num_options, argc, argv);
  }
    break;
*/
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }

  // Get required arguments
  const char* protein_database_name = get_string_parameter_pointer("protein-database");
  const char* target_input = get_string_parameter_pointer("target input");
  // Prepare the output files.
  OutputFiles output(application);
  
  // Perform the analysis.
  MatchCollection* match_collection = NULL;
  switch(command) {
  case QVALUE_COMMAND:
    match_collection = run_qvalue(target_input,
      protein_database_name,
                                  output);
    break;
/*
  case PERCOLATOR_COMMAND:
    match_collection = run_percolator(input_directory,
                                      protein_database_name,
                                      output);
    break;
*/
  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;
  }
  delete match_collection;
  output.writeFooters();
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  string name = application->getName();

  carp(CARP_INFO, "Finished crux %s.", name.c_str());
  delete application;
}

#ifdef _MSC_VER
// The Microsoft 10.0 C++ compiler has trouble resolving the proper virtual
// function call when the STL make_pair is combined with the STL ptr_fun.
// They promise to fix this in v11, but until then we create our own wrapper
// for this use of make_pair. (See corresponding ifdef block in compute_PEP)
pair<double,bool> make_pair(double db, bool b) {
    return std::pair<double,bool>(db, b);
}
#endif

/**
 * Compute posterior error probabilities (PEP) from the given target
 * and decoy scores.
 * \returns A newly allocated array of PEP for the target scores
 * sorted.
 */
double* compute_PEP(double* target_scores, ///< scores for target matches
                    int num_targets,       ///< size of target_scores
                    double* decoy_scores,  ///< scores for decoy matches
                    int num_decoys,         ///< size of decoy_scores
                    bool ascending ///< are the scores ascending or descending
){
  if( target_scores == NULL || decoy_scores == NULL 
      || num_targets == 0 || num_decoys == 0 ){
    carp(CARP_FATAL, "Cannot compute PEP without target or decoy scores.");
  }

  // put all of the scores in a single vector of pairs: score, is_target
  vector<pair<double, bool> > score_labels;
#ifdef _MSC_VER
  // There is a bug in Microsoft's implementation of
  // make_pair<> that keeps this code from working.
  // They promise to fix it in VC 11
  // https://connect.microsoft.com/VisualStudio/feedback/details/606746/incorrect-overload-resolution
  // FIXME: find workaround until then
  transform(target_scores, target_scores + num_targets,
            back_inserter(score_labels),
            bind2nd(ptr_fun<double,bool,pair<double, bool> >(make_pair), true));
  transform(decoy_scores, decoy_scores + num_decoys,
            back_inserter(score_labels),
            bind2nd(ptr_fun<double,bool,pair<double, bool> >(make_pair), false));

#else
  transform(target_scores, target_scores + num_targets,
            back_inserter(score_labels),
            bind2nd(ptr_fun(make_pair<double, bool>), true));
  transform(decoy_scores, decoy_scores + num_decoys,
            back_inserter(score_labels),
            bind2nd(ptr_fun(make_pair<double, bool>), false));
#endif

  // sort them 
  if (ascending) {
    sort(score_labels.begin(), score_labels.end());
    PosteriorEstimator::setReversed(true);
  } else {
    sort(score_labels.begin(), score_labels.end(),
       greater<pair<double, bool> > ());  
  }
  // get p-values
  vector<double> pvals;
  PosteriorEstimator::getPValues(score_labels, pvals);
  
  // estimate pi0
  double pi0 = PosteriorEstimator::estimatePi0(pvals);

  // estimate PEPs
  vector<double> PEP_vector;
  PosteriorEstimator::estimatePEP(score_labels, pi0, PEP_vector, 
                                  true);  // include decoy PEPs

  // now score_labels and PEPs are similarly sorted

  // pull out the PEPs in the order that the scores were given
  double* PEP_array = new double[PEP_vector.size()];

  for(int target_idx = 0; target_idx < num_targets; target_idx++){
    
    // the score to return next    
    double curr_target_score = target_scores[target_idx];

    // find its position in score_labels
    vector< pair<double, bool> >::iterator found_score_pos;
    if (ascending) {
      found_score_pos 
        = lower_bound(score_labels.begin(), score_labels.end(), 
                      make_pair(curr_target_score, true));
    } else {
      found_score_pos 
        = lower_bound(score_labels.begin(), score_labels.end(), 
                    make_pair(curr_target_score, true),
                    greater<pair<double, bool> >()); 
    }

    size_t found_index = distance(score_labels.begin(), found_score_pos);

    // pull out the PEP at the same position in PEP_vector
    PEP_array[target_idx] = PEP_vector[found_index];
   }

  return PEP_array;
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
MatchCollection* run_percolator(
  char* input_directory, 
  char* fasta_file, // actually fasta or index
  OutputFiles& output){ 
  input_directory= NULL; 
  fasta_file=NULL; 
  output= NULL; 

/*
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

  int num_target_matches = 0;
  int num_decoy_matches = 0;  // in first decoy set
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
      num_target_matches = match_collection->getMatchTotal();
      
      // Call that initiates percolator.
      pcInitiate(
                 (NSet)match_collection_iterator->getNumberCollections(), 
                 NUM_FEATURES, 
                 num_spectra, 
                 feature_names, 
                 pi0);
      
      free(num_spectra);

      // Call that sets verbosity level
      // 0 is quiet, 2 is default, 5 is more than you want
      if(verbosity < CARP_ERROR){
        pcSetVerbosity(0);
      }    
      else if(verbosity < CARP_INFO){
        pcSetVerbosity(1); // FIXME
      }
      else{
        pcSetVerbosity(5); // FIXME
      }
      
    } else if( set_idx == 1 ){
      num_decoy_matches = match_collection->getMatchTotal();
    }

    // create iterator, to register each PSM feature.
    carp(CARP_DETAILED_DEBUG, "Registering PSMs");
    match_iterator = new MatchIterator(match_collection, XCORR, false);
    
    while(match_iterator->hasNext()){
      match = match_iterator->next();
      // Register PSM with features
      features = match->getPercolatorFeatures(match_collection);

      output.writeMatchFeatures(match, features, NUM_FEATURES);
      pcRegisterPSM((SetType)set_idx, 
                    NULL, // no sequence used
                    features);
      
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

  double* decoy_scores = new double[num_decoy_matches]; //first decoy set
  double* PEPs = NULL;

  carp(CARP_DETAILED_DEBUG, "Processing PSMs targets:%i decoys:%i", num_target_matches, num_decoy_matches);
  // Start processing
  pcExecute(); 
  pcGetScores(results_score, results_q); 
  pcGetDecoyScores(decoy_scores);
  PEPs = compute_PEP(results_score, num_target_matches,
                     decoy_scores, num_decoy_matches);
  target_match_collection->fillResult(results_q, PERCOLATOR_QVALUE, true);
  target_match_collection->fillResult(PEPs, PERCOLATOR_PEP, true);
  target_match_collection->fillResult(results_score, PERCOLATOR_SCORE, false);
  pcCleanUp();
  
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
  delete[] decoy_scores;
  delete[] PEPs;

  // TODO put free back in. took out because glibc claimed it was corrupted
  // double linked list
  // free_parameters();

  return target_match_collection;
*/
return NULL;
}





