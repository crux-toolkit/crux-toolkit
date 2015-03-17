/**
 * \file analyze_psms.cpp
 */

#include <vector>
#include <functional>
#ifdef _MSC_VER
#include <iterator>
#endif
#include "analyze_psms.h"
#include "app/ComputeQValues.h"
#include "app/qranker-barista/QRanker.h"
#include "model/MatchCollectionIterator.h"
#include "io/MatchCollectionParser.h"
#include "model/MatchIterator.h"
#include "PosteriorEstimator.h"
#include "util/Params.h"

/**
 * \brief Takes a directory containing PSM files and a protein index
 * and analyzes the PSMs using compute-q-values, percolator or q-ranker.
 */
void analyze_matches_main(
  int argc,
  char** argv
){

  CruxApplication* application = NULL;

  application = new ComputeQValues();
  application->initialize(argc, argv);

  // Get required arguments
  vector<string> input_files = Params::GetStrings("target input");  

  if (get_boolean_parameter("list-of-files")){
    get_files_from_list(input_files[0], input_files);
  }
  
  // Prepare the output files.
  OutputFiles output(application);
  COMMAND_T command;
  if (get_string_parameter("estimation-method") == "tdc") {
    command = TDC_COMMAND;
  } else {
    command = MIXMAX_COMMAND;
  }
  
  // Perform the analysis.
  MatchCollection* match_collection = run_qvalue(input_files, output, command);

  delete match_collection;
  output.writeFooters();

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
//  pi0 = estimate_pi0(target_scores, num_targets, decoy_scores, num_decoys, ascending);

  // put all of the scores in a single vector of pairs: score, is_target
  vector<pair<double, bool> > score_labels;

  transform(target_scores, target_scores + num_targets,
            back_inserter(score_labels),
            bind2nd(ptr_fun<double,bool,pair<double, bool> >(make_pair), true));
  transform(decoy_scores, decoy_scores + num_decoys,
            back_inserter(score_labels),
            bind2nd(ptr_fun<double,bool,pair<double, bool> >(make_pair), false));

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

return NULL;
}





