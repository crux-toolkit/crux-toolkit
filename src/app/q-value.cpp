/*************************************************************************//**
 * \file q-value.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * \brief  Given as input a directory containing binary psm files,
 * a protein database, and an optional parameter file, analyze the
 * matches (with percolator or q-value) and return scores indicating
 * how good the matches are. 
 *
 * Handles at most 4 files (target and decoy).  Expects psm files to
 * start with <fileroot>.se and 
 * end with the extension '.txt' and decoys to end with
 * '-decoy#.txt'.  Multiple target files in the given directory are
 * concatinated together and presumed to be non-overlaping parts of
 * the same ms2 file. 
 ****************************************************************************/
#include "q-value.h"
#include "io/MatchCollectionParser.h"
#include "analyze_psms.h"
#include "PosteriorEstimator.h"
#include "util/crux-file-utils.h"

#include <map>
#include <utility>

using namespace std;
using namespace Crux;

static const int MAX_PSMS = 10000000;
// 14th decimal place
static const double EPSILON = 0.00000000000001;

/**
* Find the best-scoring match for each peptide in a given collection.
* Only consider the top-ranked PSM per spectrum.
*
* Results are stored in the given match collection.
*/
static void identify_best_psm_per_peptide
(MatchCollection* all_matches,
 SCORER_TYPE_T score_type)
{
  /* Instantiate a hash table.  key = peptide; value = maximal xcorr
     for that peptide. */
  map<string, FLOAT_T> best_score_per_peptide;

  // Store in the hash the best score per peptide.
  MatchIterator* match_iterator 
    = new MatchIterator(all_matches, score_type, false);
  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();

    // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char *peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position == best_score_per_peptide.end()) {
        best_score_per_peptide[peptide] = this_score;
      } else {
        // FIXME: Need a generic compare operator for score_type.
        if (map_position->second < this_score) {
          best_score_per_peptide[peptide] = this_score;
        }
      }
      free(peptide);
    }
  }
  delete match_iterator;


  // Set the best_per_peptide Boolean in the match, based on the hash.
  match_iterator = new MatchIterator(all_matches, score_type, false);
  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();

     // Skip matches that are not top-ranked.
    if (match->getRank(score_type) == 1) {
      char* peptide = match->getModSequenceStrWithSymbols();
      FLOAT_T this_score = match->getScore(score_type);

      map<string, FLOAT_T>::iterator map_position 
        = best_score_per_peptide.find(peptide);

      if (map_position->second == this_score) {
        match->setBestPerPeptide();
        
        // Prevent ties from causing two peptides to be best.
        best_score_per_peptide[peptide] = HUGE_VAL;
      }
      
      free(peptide);
    }
  }
  delete match_iterator;
}


/**
 * The q-value is defined as the minimum FDR at which a given score is
 * deemed significant.  This function takes a list of FDRs and
 * converts them into q-values.  The FDRs should be ordered from
 * lowest to highest, sorted according to the underlying score.
 */
static void convert_fdr_to_qvalue 
  (FLOAT_T* qvalues,     ///< Come in as FDRs, go out as q-values.
   int      num_values)
{
  FLOAT_T prev_fdr = qvalues[num_values - 1];
  int idx;
  for (idx=num_values - 2; idx >= 0; idx--){
    carp(CARP_DETAILED_DEBUG, "fdr[%i] = %.10f", idx, qvalues[idx]);
    FLOAT_T this_fdr = qvalues[idx];
    if (prev_fdr < this_fdr) {
      qvalues[idx] = prev_fdr;
    }
    prev_fdr = qvalues[idx];
    carp(CARP_DETAILED_DEBUG, "qvalue[%i] = %.10f", idx, qvalues[idx]);
  }
}

/**
 * Store two parallel arrays of floats in a hash table.
 *
 * The new hash table must be freed by the caller.
 */
map<FLOAT_T, FLOAT_T>* store_arrays_as_hash
  (FLOAT_T* keys, 
   FLOAT_T* values,
   int      num_values
){

  map<FLOAT_T, FLOAT_T>* return_value = new map<FLOAT_T, FLOAT_T>();

  int idx;
  for (idx=0; idx < num_values; idx++){
    carp(CARP_DETAILED_DEBUG, "%g maps to %g", keys[idx], values[idx]);
    (*return_value)[keys[idx]] = values[idx];
  }
  return(return_value);
}

/**
 * \brief Compute q-values from a given set of scores, using a second
 * set of scores as an empirical null.  Sorts the incoming target
 * scores and returns a corresponding list of q-values.
 *
 * This function is only exported to allow unit testing.
 */
FLOAT_T* compute_decoy_qvalues_tdc(
  FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys,
  bool     forward,
  FLOAT_T  pi_zero
){
  if ((num_targets == 0) || (num_decoys == 0)) {
    carp(CARP_FATAL, "Cannot compute q-values (%d targets, %d nulls).",
         num_targets, num_decoys);
  }
  carp(CARP_DEBUG, "Computing decoy q-values.");

  int target_idx;
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    carp(CARP_DEBUG, "target_scores[%d]=%g decoy_scores[%d]=%g",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  // Sort both sets of scores.
  if (forward) {
    sort(target_scores, target_scores + num_targets);
    sort(decoy_scores, decoy_scores + num_decoys);    
  } else {
    sort(target_scores, target_scores + num_targets, greater<FLOAT_T>());
    sort(decoy_scores, decoy_scores + num_decoys, greater<FLOAT_T>());
  }
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    carp(CARP_DEBUG, "target_scores[%d]=%g decoy_scores[%d]=%g",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  // Precompute the ratio of targets to decoys.
  FLOAT_T targets_to_decoys = (FLOAT_T)num_targets / (FLOAT_T)num_decoys;

  // Compute false discovery rate for each target score.
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_targets, sizeof(FLOAT_T));
  int decoy_idx = 0;
  for (target_idx = 0; target_idx < num_targets; target_idx++) {
    FLOAT_T target_score = target_scores[target_idx];

    // Find the index of the first decoy score greater than this target score.
    if (forward) {
      while ((decoy_idx < num_decoys) &&
             (decoy_scores[decoy_idx] < target_score)) {
        decoy_idx++;
      }
    } else {   
      while ((decoy_idx < num_decoys) &&
             (decoy_scores[decoy_idx] > target_score)) {
        decoy_idx++;
      }
    }

    // FDR = #decoys / #targets
    FLOAT_T fdr = /*pi_zero * targets_to_decoys * */
      ((FLOAT_T)decoy_idx / (FLOAT_T)(target_idx + 1));
    carp(CARP_DEBUG, "target_idx=%d target_score=%g decoy_idx=%d fdr=%g",
         target_idx, target_score, decoy_idx, fdr);
    
    if ( fdr > 1.0 ){
      fdr = 1.0;
    }
    
    qvalues[target_idx] = fdr;
  }
  
  // Convert the FDRs into q-values.
  convert_fdr_to_qvalue(qvalues, num_targets);

  return (qvalues);
}

FLOAT_T estimate_pi0( FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys,
  bool     ascending
){
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
  
  // estimate pi_zero
  FLOAT_T pi_zero = PosteriorEstimator::estimatePi0(pvals);

  carp(CARP_INFO, "Estimated pi_zero = %f", pi_zero);
  return pi_zero;
}

//#ifdef _MSC_VER
// The Microsoft 10.0 C++ compiler has trouble resolving the proper virtual
// function call when the STL make_pair is combined with the STL ptr_fun.
// They promise to fix this in v11, but until then we create our own wrapper
// for this use of make_pair. (See corresponding ifdef block in compute_PEP)
pair<double,bool> make_pair(double db, bool b) {
    return std::pair<double,bool>(db, b);
}
//#endif

/**
 * \brief Compute q-values using mix-max procedure. This part is a
 * reimplementation of Uri Keich's code written in R.
 *
 */
FLOAT_T* compute_decoy_qvalues_mixmax(
  FLOAT_T* target_scores,
  int      num_targets,
  FLOAT_T* decoy_scores,
  int      num_decoys,
  bool     ascending,
  FLOAT_T  pi_zero
){
  if ((num_targets == 0) || (num_decoys == 0)) {
    carp(CARP_FATAL, "Cannot compute q-values (%d targets, %d decoys).",
         num_targets, num_decoys);
  }
  if (num_targets != num_decoys) {
    carp(CARP_FATAL, "Mix-Max requires equal number of decoy and target scores "
                     "(%d targets, %d decoys) from separate target-decoy search.",
         num_targets, num_decoys);
  }
  //estimate pi0 from data if it is not given.
  if (pi_zero == 1.0) {
//    pi_zeros = estimate_pi0(target_scores, num_targets, decoy_scores, num_decoys, ascending);
    
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
      
      // estimate pi_zero
      pi_zero = PosteriorEstimator::estimatePi0(pvals);

      carp(CARP_INFO, "Estimated pi_zero = %f", pi_zero);
 //     pi_zero =0.5302161;

  }
  // continue with mix-max procedure
  int target_idx;
  for (target_idx = 0; target_idx < num_targets; ++target_idx) {
    carp(CARP_DEBUG, "target_scores[%d]=%lf decoy_scores[%d]=%lf",
         target_idx, target_scores[target_idx],
         target_idx, decoy_scores[target_idx]);
  }

  //Sort decoy and target stores
  ascending = false;
  if (ascending) {
    sort(target_scores, target_scores + num_targets);
    sort(decoy_scores, decoy_scores + num_decoys);    
  } else {
    sort(target_scores, target_scores + num_targets, greater<FLOAT_T>());
    sort(decoy_scores, decoy_scores + num_decoys, greater<FLOAT_T>());
  }

  //histogram of the target scores.

  vector<double> z_hist;  
  z_hist.reserve(num_decoys+1);
  int idx = 0;
  int cnt;
  int i;
  for (i = 0; i < num_decoys; ++i) {
    cnt = 0;
    while (ascending ? 
        target_scores[idx] < decoy_scores[i] : 
        target_scores[idx] >= decoy_scores[i]) {
      ++cnt;
      ++idx;
    }
    z_hist[i] = (double)cnt;///(double)num_targets;
  }
  z_hist[num_decoys] = (double)(num_targets - idx);///num_targets;

  for (i = 1; i <= num_decoys; ++i){
    z_hist[i] += z_hist[i-1];
  }

  vector<double> p_T_and_X_lt_Y;
  p_T_and_X_lt_Y.reserve(num_targets);
  double estPx_lt_zj;

  for (i = 0; i < num_targets; ++i){
    estPx_lt_zj = (double)( z_hist[i+1] - pi_zero * ((double)i+0.5) ) / (double)(1.0-pi_zero) / (i+0.5);
    estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
    estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
    p_T_and_X_lt_Y[i] = estPx_lt_zj * ((1.0-pi_zero));    
  }
  
  FLOAT_T* fdrmod = new FLOAT_T[num_targets];
  double E_f1_mod_run_tot = 0;
  double n_z_gt_w = 0;
  int j = num_targets-1;
  double qvalue;
  for (i = num_targets-1; i >= 0; --i){ 

    while (j >= 0 && (ascending ? decoy_scores[j] < target_scores[i] : decoy_scores[j] < target_scores[i])) { 
      E_f1_mod_run_tot  += p_T_and_X_lt_Y[j];
      ++n_z_gt_w;
      --j;
    }
    //when j < 0 the running sums remain frozen till the end of the loop
    qvalue = (n_z_gt_w * pi_zero + E_f1_mod_run_tot) / (num_targets - i + 0.5-1);
    fdrmod[i] = qvalue > 1 ? 1 : qvalue;  
  }
  return fdrmod;
}


/**
 * Use the given p-values to estimate PEP using the
 * PosteriorEstimator.
 * \returns A newly allocated array of PEP values.
 */
FLOAT_T* compute_PEP_from_pvalues(FLOAT_T* pvalues, int num_pvals){

  // convert the -ln(pval) to pval
  vector<double> pvalues_vector(pvalues, pvalues + num_pvals);
  for(size_t val_idx = 0; val_idx < pvalues_vector.size(); val_idx++){
    pvalues_vector[val_idx] = exp(-pvalues_vector[val_idx]);
  }

  // sort them
  sort(pvalues_vector.begin(), pvalues_vector.end());

  // put them in a score_label vector
  vector<pair<double, bool> > score_label;

#ifdef _MSC_VER
  // There is a bug in Microsoft's implementation of
  // make_pair<> that keeps this code from working.
  // They promise to fix it in VC 11
  // https://connect.microsoft.com/VisualStudio/feedback/details/606746/incorrect-overload-resolution
  score_label.reserve(pvalues_vector.size());
  for (vector<double>::const_iterator i = pvalues_vector.begin();
       i != pvalues_vector.end();
       i++) {
    score_label.push_back(make_pair(*i, true));
  }
#else
  transform(pvalues_vector.begin(),
            pvalues_vector.end(),
            back_inserter(score_label),
            bind2nd(ptr_fun(make_pair<double,bool>), true));
#endif

  // create decoy p-values
  double step = 1.0 / 2.0 / (double)num_pvals;
  for(int val_idx = 0; val_idx < num_pvals; val_idx++){
    score_label.push_back(make_pair<double, bool>(step * (1 + 2 * val_idx),
                                               false));
  }

  // sort ascending order
  sort(score_label.begin(), score_label.end());
  PosteriorEstimator::setReversed(true);

  // estimate PEPs 
  double pi0 = PosteriorEstimator::estimatePi0(pvalues_vector);
  vector<double> PEP_vector;
  PosteriorEstimator::estimatePEP(score_label, pi0, PEP_vector );

  // return values
  FLOAT_T* PEPs = new FLOAT_T[PEP_vector.size()];
  for(size_t pep_idx = 0; pep_idx < PEP_vector.size(); pep_idx++){
    PEPs[pep_idx] = PEP_vector[pep_idx];
    carp(CARP_DEBUG, "pep[%i]=%f", pep_idx, PEPs[pep_idx]);
  }
  return PEPs;
}

/**
 * \brief Compute a q-values based on what is in the PSM files in the
 * directory.  Store q-values in the match collection returned.
 *
 * If p-values were computed, then perform Benjamini-Hochberg q-value
 * calculations. Otherwise, if decoys are present, then rank on xcorr
 * and compute empirical q-values based on the number of decoys and
 * targets above the score threshold.
 *
 * \returns a collection of target PSMs with one q-value and one PEP
 * in each match.
 */
MatchCollection* run_qvalue(
  vector<string>& input_files,
  const string& fasta_file,
  OutputFiles& output,
  COMMAND_T command  
  ){
  
  if (input_files.size() == 0) {
    carp(CARP_FATAL, "No search paths found!");
  }
  
  bool ascending = get_boolean_parameter("smaller-is-better");
  SCORER_TYPE_T score_type = INVALID_SCORER_TYPE;
  
  if (get_string_parameter("score") == "exact p-value") {
    score_type = TIDE_SEARCH_EXACT_PVAL;
  } else if (get_string_parameter("score") == "xcorr score") {        
    score_type = XCORR;
  }

  // Create two match collections, for targets and decoys.
  MatchCollection* decoy_matches = new MatchCollection();
  MatchCollection* target_matches = new MatchCollection();

  bool distinct_matches = false; 
//  for (vector<string>::iterator iter = input_files.begin(); iter != input_files.end(); ++iter) {

vector<string>::iterator iter = input_files.begin();
    string target_path = *iter;
    string decoy_path = *iter;
  
    check_target_decoy_files(target_path, decoy_path);

    if (!file_exists(target_path)) {
      carp(CARP_FATAL, "Target file %s not found", target_path.c_str());
    }
    
    if (!file_exists(decoy_path)) {
      if (command == MIXMAX_COMMAND) {
        carp(CARP_FATAL, "Decoy file separate target-decoy search is required for q-value calculation");
      }
      carp(CARP_DEBUG, "Decoy file %s not found", decoy_path.c_str());
      decoy_path = "";    
    }

    MatchCollectionParser parser;
    MatchCollection* match_collection =
      parser.create(target_path, get_string_parameter("protein-database"));
    distinct_matches  = match_collection->getHasDistinctMatches();

    target_matches->setScoredType(TIDE_SEARCH_EXACT_PVAL,match_collection->getScoredType(TIDE_SEARCH_EXACT_PVAL));
    target_matches->setScoredType(EVALUE,match_collection->getScoredType(EVALUE));
    
    target_matches->setScoredType(DELTA_CN,match_collection->getScoredType(DELTA_CN));
    target_matches->setScoredType(SP,match_collection->getScoredType(SP));
    target_matches->setScoredType(BY_IONS_MATCHED,match_collection->getScoredType(BY_IONS_MATCHED));
    target_matches->setScoredType(BY_IONS_TOTAL,match_collection->getScoredType(BY_IONS_TOTAL));
    
    if (decoy_path != "") {
      MatchCollection* temp_collection = parser.create(decoy_path, get_string_parameter("protein-database"));
         // Mark decoy matches
      MatchIterator* temp_iter = new MatchIterator(temp_collection);
      while (temp_iter->hasNext()) {
        Crux::Match* decoy_match = temp_iter->next();
          decoy_match->setNullPeptide(true);
          if (command != TDC_COMMAND) {
            decoy_matches->addMatch(decoy_match);
          }
      }
      delete temp_iter;
      //carry out concatenated search.
      if (command == TDC_COMMAND) {
        MatchCollection* tdc_collection = new MatchCollection();
        tdc_collection->setScoredType(score_type, true);      
        MatchIterator* target_iter = new MatchIterator(match_collection);
        MatchIterator* decoy_iter = new MatchIterator(temp_collection);
        while (target_iter->hasNext() ) {
          Crux::Match* target_match = target_iter->next();
          Crux::Match* decoy_match  = decoy_iter->next();
          decoy_match->setNullPeptide(true);          
          if( target_match->getRank(XCORR) != 1 || decoy_match->getRank(XCORR) != 1 ) continue;
          if (ascending) { 
            tdc_collection->addMatch(target_match->getScore(score_type) < decoy_match->getScore(score_type) ? target_match : decoy_match);
          } else {
            tdc_collection->addMatch(target_match->getScore(score_type) > decoy_match->getScore(score_type) ? target_match : decoy_match);
          }
        }  

        delete target_iter;
        delete decoy_iter;
        delete match_collection;
        match_collection = tdc_collection;
      }
      delete temp_collection;
    }
     
    // Iterate, gathering matches into one or two collections.
    MatchIterator* match_iterator =
      new MatchIterator(match_collection, score_type, false);

    while(match_iterator->hasNext()){
      Match* match = match_iterator->next();
        
      // Only use top-ranked matches.
      if( match->getRank(XCORR) != 1 ){
        continue;
      }
      if (match->getNullPeptide() == true) {
        decoy_matches->addMatch(match);
      } else {
        target_matches->addMatch(match);
      }
      Match::freeMatch(match);
    }
    delete match_iterator;   
    delete match_collection;   
//    break;
  
  bool have_pvalues = target_matches->getScoredType(TIDE_SEARCH_EXACT_PVAL);
  bool have_evalues = target_matches->getScoredType(EVALUE);
  target_matches->setScoredType(score_type, true);
  decoy_matches->setScoredType(score_type, true);   


  // get from the input files which columns to print in the output files
  vector<bool> cols_to_print(NUMBER_MATCH_COLUMNS);
  cols_to_print[FILE_COL] = true;
  cols_to_print[SCAN_COL] = true;
  cols_to_print[CHARGE_COL] = true;
  cols_to_print[SPECTRUM_PRECURSOR_MZ_COL] = true;
  cols_to_print[SPECTRUM_NEUTRAL_MASS_COL] = true;
  cols_to_print[PEPTIDE_MASS_COL] = true;
  cols_to_print[DELTA_CN_COL] = target_matches->getScoredType(DELTA_CN);
  cols_to_print[SP_SCORE_COL] = target_matches->getScoredType(SP);
  cols_to_print[SP_RANK_COL] = target_matches->getScoredType(SP);
  cols_to_print[XCORR_SCORE_COL] = !have_pvalues;
  cols_to_print[XCORR_RANK_COL] = true;
  cols_to_print[EVALUE_COL] = have_evalues;
  cols_to_print[EXACT_PVALUE_COL] = have_pvalues;
  if (have_pvalues) {
    cols_to_print[REFACTORED_SCORE_COL] = true;
  }
  cols_to_print[BY_IONS_MATCHED_COL] = target_matches->getScoredType(BY_IONS_MATCHED);
  cols_to_print[BY_IONS_TOTAL_COL] = target_matches->getScoredType(BY_IONS_TOTAL);

  if (distinct_matches) {
    cols_to_print[DISTINCT_MATCHES_SPECTRUM_COL] = true;
  } else {
    cols_to_print[MATCHES_SPECTRUM_COL] = true;
  }

  cols_to_print[SEQUENCE_COL] = true;
  cols_to_print[CLEAVAGE_TYPE_COL] = true;
  cols_to_print[PROTEIN_ID_COL] = true;
  cols_to_print[FLANKING_AA_COL] = true;
  cols_to_print[ETA_COL] = have_pvalues;
  cols_to_print[BETA_COL] = have_pvalues;
  cols_to_print[SHIFT_COL] = have_pvalues;
  cols_to_print[CORR_COL] = have_pvalues;

  output.writeHeaders(cols_to_print);

  // Compute q-values from p-values.
  FLOAT_T* pvalues = NULL; // N.B. Misnamed for decoy calculation.
  int num_pvals = target_matches->getMatchTotal();
  FLOAT_T* qvalues = NULL;
  FLOAT_T* decoy_scores;
  int num_decoys  ;
  num_decoys = decoy_matches->getMatchTotal();
  carp(CARP_INFO,
       "There are %d target and %d decoy PSMs for q-value computation.",
       num_pvals, num_decoys);
  decoy_scores = NULL;

  pvalues = target_matches->extractScores(score_type);
  decoy_scores = decoy_matches->extractScores(score_type);
  switch (command){
    case TDC_COMMAND:

      qvalues = compute_decoy_qvalues_tdc(pvalues, num_pvals, 
                                      decoy_scores, num_decoys, 
                                      get_boolean_parameter("smaller-is-better"),
                                      1.0);

      break;
    case MIXMAX_COMMAND:
 
      qvalues = compute_decoy_qvalues_mixmax(pvalues, num_pvals, 
                                      decoy_scores, num_decoys,
                                      get_boolean_parameter("smaller-is-better"), 
                                      get_double_parameter("pi-zero"));
    
      break;
   }
   free(decoy_scores);
    
  // Store p-values to q-values as a hash, and then assign them.
  map<FLOAT_T, FLOAT_T>* qvalue_hash 
    = store_arrays_as_hash(pvalues, qvalues, num_pvals);

  target_matches->assignQValues(qvalue_hash, score_type);

  free(pvalues);
  free(qvalues);
  delete qvalue_hash;
  // Identify PSMs that are top-scoring per peptide.
  identify_best_psm_per_peptide(target_matches, score_type);
  
  // Store targets by score.
  target_matches->sort(score_type);
  output.writeMatches(target_matches);

  delete decoy_matches;

  return(target_matches);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
