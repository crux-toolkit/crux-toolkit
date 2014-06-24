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
#include "MatchCollectionParser.h"
#include "analyze_psms.h"
#include "PosteriorEstimator.h"

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
 * Use the Benjamini-Hochberg procedure to convert a given set of
 * p-values into q-values.  
 *
 * Assumes that the input is an array of negative log p-values.  The
 * output q-values are not log-transformed.
 *
 * This function uses the command line parameter "pi-zero".
 */
FLOAT_T* compute_qvalues_from_pvalues(
  FLOAT_T* pvalues, 
  int      num_pvals,
  FLOAT_T  pi_zero
){

  // sort the - log p-values in descending order
  sort(pvalues, pvalues + num_pvals, greater<FLOAT_T>());

  // convert the p-values into FDRs using Benjamini-Hochberg
  FLOAT_T* qvalues = (FLOAT_T*)mycalloc(num_pvals, sizeof(FLOAT_T));
  int idx;
  for (idx=0; idx < num_pvals; idx++){
    carp(CARP_DETAILED_DEBUG, "pvalue[%i] = %.10f", idx, exp(-pvalues[idx]));
    double fdr = (exp(-pvalues[idx]) / (idx + 1)) 
      * (FLOAT_T)num_pvals * pi_zero;
    qvalues[idx] = fdr;
    carp(CARP_DETAILED_DEBUG, "FDR[%i] = %.10f", idx, qvalues[idx]);
  }

  // convert the FDRs into q-values
  convert_fdr_to_qvalue(qvalues, num_pvals);

  return(qvalues);
}

/**
 * \brief Compute q-values from a given set of scores, using a second
 * set of scores as an empirical null.  Sorts the incoming target
 * scores and returns a corresponding list of q-values.
 *
 * This function is only exported to allow unit testing.
 */
FLOAT_T* compute_decoy_qvalues(
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
    FLOAT_T fdr = pi_zero * targets_to_decoys * 
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

  return(qvalues);
}

/**
 * A wrapper to take care of the FLOAT_T to double conversion.  Calls
 * compute_PEP in analyse_psms.
 * \returns A newly allocated array of PEP values sorted in the same
 * order as the targets.
 */
FLOAT_T* compute_PEP_local(FLOAT_T* targets,
                          int num_targets, 
                          FLOAT_T* decoys, 
                          int num_decoys,
                          bool forward){

  double* targets_d = new double[num_targets];
  for(int val_idx = 0; val_idx < num_targets; val_idx++){
    targets_d[val_idx] = targets[val_idx];
  }

  double* decoys_d = new double[num_decoys];
  for(int val_idx = 0; val_idx < num_decoys; val_idx++){
    decoys_d[val_idx] = decoys[val_idx];
  }

  double* PEPs_d = compute_PEP(targets_d, num_targets, decoys_d, num_decoys);

  FLOAT_T* PEPs = new FLOAT_T[num_targets];
  for(int val_idx = 0; val_idx < num_targets; val_idx++){
    PEPs[val_idx] = PEPs_d[val_idx];
  }

  delete []targets_d;
  delete []decoys_d;
  delete []PEPs_d;
  return PEPs;
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
  const char* input_file, 
  const char* fasta_file,
  OutputFiles& output 
  ){


  string target_path = string(input_file);
  string decoy_path = string(input_file);
  
  check_target_decoy_files(target_path, decoy_path);

  if (!file_exists(target_path)) {
    carp(CARP_FATAL, "Target file %s not found", target_path.c_str());
  }
  
  if (!file_exists(decoy_path)) {
    carp(CARP_DEBUG, "Decoy file %s not found", decoy_path.c_str());
    decoy_path = "";
  }

  bool have_decoys = false;
  MatchCollectionParser parser;
  MatchCollection* match_collection =
    parser.create(target_path.c_str(), get_string_parameter_pointer("protein-database"));
  bool distinct_matches = match_collection->getHasDistinctMatches();

  MatchCollection* decoy_matches = new MatchCollection();
  // Create two match collections, for targets and decoys.
  MatchCollection* target_matches = new MatchCollection();
  

  if (decoy_path != "") {
    MatchCollection* temp_collection = parser.create(decoy_path.c_str(), get_string_parameter_pointer("protein-database"));
       // Mark decoy matches
    MatchIterator* temp_iter = new MatchIterator(temp_collection);
    while (temp_iter->hasNext()) {
      Crux::Match* decoy_match = temp_iter->next();
      if (decoy_match->getRank(XCORR) == 1) {
        decoy_match->setNullPeptide(true);
        decoy_matches->addMatch(decoy_match);
        have_decoys = true;
      }
    }
    delete temp_iter;
    delete temp_collection;
  }
  
  target_matches->setScoredType(XCORR, true);
  decoy_matches->setScoredType(XCORR, true);

  // Did we find something from which to get q-values?
  bool have_pvalues = match_collection->getScoredType(LOGP_BONF_WEIBULL_XCORR);
  bool have_evalues = match_collection->getScoredType(EVALUE);
  
  target_matches->setScoredType(EVALUE, have_evalues);
  
  // Iterate, gathering matches into one or two collections.
  MatchIterator* match_iterator =
    new MatchIterator(match_collection, XCORR, false);
  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
      
    // Only use top-ranked matches.
    if( match->getRank(XCORR) != 1 ){
      continue;
    }
    if (match->getNullPeptide() == true) {
      decoy_matches->addMatch(match);
      have_decoys = true;
    } else {
      target_matches->addMatch(match);
    }
    Match::freeMatch(match);
  }
  delete match_iterator;

  // get from the input files which columns to print in the output files
  vector<bool> cols_to_print(NUMBER_MATCH_COLUMNS);
  cols_to_print[FILE_COL] = true;
  cols_to_print[SCAN_COL] = true;
  cols_to_print[CHARGE_COL] = true;
  cols_to_print[SPECTRUM_PRECURSOR_MZ_COL] = true;
  cols_to_print[SPECTRUM_NEUTRAL_MASS_COL] = true;
  cols_to_print[PEPTIDE_MASS_COL] = true;
  cols_to_print[DELTA_CN_COL] = match_collection->getScoredType(DELTA_CN);
  cols_to_print[SP_SCORE_COL] = match_collection->getScoredType(SP);
  cols_to_print[SP_RANK_COL] = match_collection->getScoredType(SP);
  cols_to_print[XCORR_SCORE_COL] = true;
  cols_to_print[XCORR_RANK_COL] = true;
  cols_to_print[EVALUE_COL] = have_evalues;
  cols_to_print[PVALUE_COL] = have_pvalues;
  cols_to_print[BY_IONS_MATCHED_COL] = match_collection->getScoredType(BY_IONS_MATCHED);
  cols_to_print[BY_IONS_TOTAL_COL] = match_collection->getScoredType(BY_IONS_TOTAL);

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
  delete match_collection;


  output.writeHeaders(cols_to_print);

  // Compute q-values from p-values.
  FLOAT_T* pvalues = NULL; // N.B. Misnamed for decoy calculation.
  int num_pvals = target_matches->getMatchTotal();
  FLOAT_T* qvalues = NULL;
  FLOAT_T* PEPs = NULL;
  SCORER_TYPE_T score_type = INVALID_SCORER_TYPE;
  if (have_pvalues == true) {
    carp(CARP_DEBUG, "There are %d PSMs for q-value computation.", num_pvals);
    target_matches->setScoredType(LOGP_BONF_WEIBULL_XCORR, 
                                     true);
    pvalues = target_matches->extractScores(LOGP_BONF_WEIBULL_XCORR);
    qvalues = compute_qvalues_from_pvalues(pvalues, num_pvals,
                                           get_double_parameter("pi-zero"));
    PEPs = compute_PEP_from_pvalues(pvalues, num_pvals);
    score_type = LOGP_BONF_WEIBULL_XCORR;
  }

  // Compute q-values from the XCorr decoy distribution.
  else if (have_decoys == true) {
    int num_decoys = decoy_matches->getMatchTotal();
    carp(CARP_INFO,
         "There are %d target and %d decoy PSMs for q-value computation.",
         num_pvals, num_decoys);
    FLOAT_T* decoy_scores = NULL;

    pvalues = target_matches->extractScores(XCORR);
    decoy_scores = decoy_matches->extractScores(XCORR);
    score_type = XCORR;

    qvalues = compute_decoy_qvalues(pvalues, num_pvals, 
                                    decoy_scores, num_decoys, false, 
                                    get_double_parameter("pi-zero"));


    PEPs = compute_PEP_local(pvalues, num_pvals, decoy_scores, num_decoys, false);

    free(decoy_scores);
  }
  // Fatal: Cannot compute q-values.
  else {
    carp(CARP_FATAL, "Cannot compute q-values without decoy PSMs or p-values.");
  }

  
  
  // Store p-values to q-values as a hash, and then assign them.
  map<FLOAT_T, FLOAT_T>* qvalue_hash 
    = store_arrays_as_hash(pvalues, qvalues, num_pvals);

  target_matches->assignQValues(qvalue_hash, score_type);

  // Store p-values to PEP as a has and then assign them
  map<FLOAT_T, FLOAT_T>* PEP_hash 
        = store_arrays_as_hash(pvalues, PEPs, num_pvals);
  target_matches->assignPEPs(PEP_hash, score_type);

  free(pvalues);
  free(qvalues);
  delete PEPs;
  delete qvalue_hash;
  delete PEP_hash;

  if (have_evalues) {
    pvalues = target_matches->extractScores(EVALUE);
    FLOAT_T* decoy_scores = decoy_matches->extractScores(EVALUE);
    int num_decoys = decoy_matches->getMatchTotal();
    
    score_type = EVALUE;
    qvalues = compute_decoy_qvalues(pvalues, num_pvals, 
                                    decoy_scores, num_decoys, true, 
                                    get_double_parameter("pi-zero"));
    PEPs = compute_PEP_local(pvalues, num_pvals, decoy_scores, num_decoys, true);
    free(decoy_scores);
    qvalue_hash = store_arrays_as_hash(pvalues, qvalues, num_pvals);
    target_matches->assignQValues(qvalue_hash, score_type);
    PEP_hash = store_arrays_as_hash(pvalues, PEPs, num_pvals);
    target_matches->assignPEPs(PEP_hash, score_type);
    
    free(pvalues);
    free(qvalues);
    delete PEPs;
    delete qvalue_hash;
    delete PEP_hash;
  }
  
  // Identify PSMs that are top-scoring per peptide.
  identify_best_psm_per_peptide(target_matches, score_type);

  // Compute peptide-level q-values.
  //  compute_decoy_q_values(all_matches/
  // true); // Do peptide-level scoring.

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
