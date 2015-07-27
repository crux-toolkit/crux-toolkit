/*************************************************************************//**
 * \file Match.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * \brief Object for matching a peptide and a spectrum, generate
 * a preliminary score(e.g., Sp) 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <set>
#include <map>
#include "io/carp.h"
#include "ProteinPeptideIterator.h"
#include "Spectrum.h"
#include "Ion.h"
#include "IonSeries.h"
#include "util/crux-utils.h"
#include "objects.h"
#include "io/OutputFiles.h"
#include "parameter.h"
#include "Scorer.h" 
#include "Match.h" 
#include "MatchCollection.h" 
#include "Peptide.h"
#include "util/Params.h"

#include <boost/filesystem.hpp>
#include <string>

#include "io/MatchFileReader.h"

using namespace std;
using namespace Crux;

vector<string> Match::file_paths_;

void Match::init() {

  spectrum_=NULL;
  peptide_=NULL;
  // initialize score, rank !!!! DEBUG
  for (unsigned int idx=0;idx<NUMBER_SCORER_TYPES;idx++) {
    match_scores_[idx] = NOT_SCORED;
    match_rank_[idx] = 0;
  }
  pointer_count_ = 0;
  b_y_ion_fraction_matched_ = 0;
  b_y_ion_matched_ = 0;
  b_y_ion_possible_ = 0;
  // default is not a null peptide match
  null_peptide_ = false;
  peptide_sequence_ = NULL;
  mod_sequence_ = NULL;
  digest_ = INVALID_DIGEST;
  post_process_match_ = 0;
  ln_experiment_size_ = 0;
  num_target_matches_ = 0;
  num_decoy_matches_ = 0;
  best_per_peptide_ = false;
  file_idx_ = -1;
}


/**
 * \returns a new memory allocated match
 */
Match::Match(){
  init();
  ++pointer_count_;
  exact_pval_search_ = false;
}

/**
 * Create a new match with the given members.
 */
Match::Match(Peptide* peptide, ///< the peptide for this match
             Spectrum* spectrum, ///< the spectrum for this match
             SpectrumZState& zstate, ///< the charge/mass of the spectrum
             bool is_decoy)///< is the peptide a decoy or not
{
  init();
  peptide_ = peptide;
  spectrum_ = spectrum;
  zstate_ = zstate;
  null_peptide_ = is_decoy;

  ++pointer_count_;
  exact_pval_search_ = false;
}

/**
 * free the memory allocated match
 * spectrum is not freed by match
 */
void Match::freeMatch(
  Match* match ///< the match to free -in
  )
{
  --match->pointer_count_;
  
  // only free match when pointer count reaches
  if(match->pointer_count_ == 0){
    delete match;
  }
}

Match::~Match() {

  // but aren't there multiple matches pointing to the same peptide?
  // if so, create a new free_shallow_match which doesn't touch the members
  if (peptide_ != NULL){
    delete peptide_;
  }
  if(post_process_match_ && spectrum_ !=NULL){
    delete spectrum_;
  }
  if (peptide_sequence_ != NULL){
    free(peptide_sequence_);
  }
  if (mod_sequence_ != NULL){
    free(mod_sequence_);
  }
}

/**
 * \brief Print the match information in sqt format to the given file
 *
 * Only crux sequest-search produces sqt files so the two scores
 * printed are always Sp and xcorr.
 */
void Match::printSqt(
  FILE* file                  ///< output stream -out
  ){

  if( file == NULL ){
    carp(CARP_ERROR, "Cannot print match to sqt file from null input");
    return;
  }

  Peptide* peptide = getPeptide();
  // this should get the sequence from the match, not the peptide
  char* sequence = getSequenceSqt();

  int b_y_total = getBYIonPossible();
  int b_y_matched = getBYIonMatched();
  
  FLOAT_T delta_cn = getScore(DELTA_CN);
  FLOAT_T score_main = getScore(XCORR);

  // write format string with variable precision
  int precision = get_int_parameter("precision");

  // print match info
  if (exact_pval_search_) {
    fprintf(file, "M\t%i\t%i\t%.*f\t%.2f\t%.*g\t%.*g\t%.*g\t%i\t%i\t%s\tU\n",
            getRank(XCORR),
            getRank(SP),
            get_int_parameter("mass-precision"),
            peptide->getPeptideMass() + MASS_PROTON,
            delta_cn,
            precision,
            getScore(TIDE_SEARCH_EXACT_PVAL),
            precision,
            getScore(TIDE_SEARCH_REFACTORED_XCORR),
            precision,
            getScore(SP),
            b_y_matched,
            b_y_total,
            sequence
            );
  } else {
    fprintf(file, "M\t%i\t%i\t%.*f\t%.2f\t%.*g\t%.*g\t%i\t%i\t%s\tU\n",
            getRank(XCORR),
            getRank(SP),
            get_int_parameter("mass-precision"),
            peptide->getPeptideMass() + MASS_PROTON,
            delta_cn,
            precision,
            score_main,
            precision,
            getScore(SP),
            b_y_matched,
            b_y_total,
            sequence
            );
  }
  free(sequence);
  
  PeptideSrc* peptide_src = NULL;
  char* protein_id = NULL;
  Protein* protein = NULL;
  
  for(PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
      iter != peptide->getPeptideSrcEnd();
      ++iter){
    string rand;
    peptide_src = *iter;
    protein = peptide_src->getParentProtein();
    protein_id = protein->getId();

    // only prepend "decoy-prefix" if we are doing a fasta search
    Database* database = protein->getDatabase();
    if( null_peptide_ 
        && (database != NULL && database->getDecoyType() == NO_DECOYS) ){
      rand = get_string_parameter("decoy-prefix"); 
    }

    // print match info (locus line), add "decoy-prefix" to locus name for decoys
    fprintf(file, "L\t%s%s\n", rand.c_str(), protein_id);      
    free(protein_id);
  }
  
  return;
}


/**
 * Print one field in the tab-delimited output file, based on column index.
 */
void Match::printOneMatchField(
  int      column_idx,             ///< Index of the column to print. -in
  MatchCollection* collection,  ///< collection holding this match -in 
  MatchFileWriter*    output_file,            ///< output stream -out
  Spectrum* spectrum,
  int      num_target_matches,            ///< num matches per spectrum -in
  int      num_decoy_matches,     ///< target matches for same spectrum -in
  int      b_y_total,              ///< total b/y ions -in
  int      b_y_matched             ///< Number of b/y ions matched. -in
) {
  switch ((MATCH_COLUMNS_T)column_idx) {
  case FILE_COL:
    if (strlen(spectrum->getFullFilename()) == 0) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, getFilePath());
    }
    else {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, spectrum->getFullFilename());
    }
    break;
  case SCAN_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     spectrum->getFirstScan());
    break;
  case CHARGE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                      getCharge());
    break;
  case SPECTRUM_PRECURSOR_MZ_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     spectrum->getPrecursorMz());
    break;
  case SPECTRUM_NEUTRAL_MASS_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getNeutralMass());
    break;
  case PEPTIDE_MASS_COL:
    {
      Peptide* peptide = getPeptide();
      double peptide_mass = peptide->getPeptideMass();
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       peptide_mass);
    }
    break;
  case DELTA_CN_COL:
    {
      FLOAT_T delta_cn = getScore(DELTA_CN);
      if( delta_cn == 0 ){// I hate -0, this prevents it
        delta_cn = 0.0;
      }
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, delta_cn);
    }
    break;
  case SP_SCORE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(SP));
    break;
  case SP_RANK_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getRank(SP));
    break;
  case XCORR_SCORE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(XCORR));
    break;
  case EXACT_PVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(TIDE_SEARCH_EXACT_PVAL));
    break;
  case REFACTORED_SCORE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(TIDE_SEARCH_REFACTORED_XCORR));
    break;
  case SIDAK_ADJUSTED_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(SIDAK_ADJUSTED));
    break;
  case XCORR_RANK_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getRank(!Params::GetBool("exact-p-value") ? XCORR : TIDE_SEARCH_EXACT_PVAL));
    break;
  case EVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getScore(EVALUE));
    break;
  case PVALUE_COL:
    {
      double log_pvalue = getScore(LOGP_BONF_WEIBULL_XCORR);
      if (P_VALUE_NA == log_pvalue) {
        output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, "NaN");
      }
      else {
        output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                         exp(-1 * log_pvalue));
      }
    }
    break;
  case WEIBULL_QVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                             getScore(LOGP_QVALUE_WEIBULL_XCORR));
    break;
  case WEIBULL_PEP_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getScore(LOGP_WEIBULL_PEP));
    break;
#ifdef NEW_COLUMNS
  case WEIBULL_PEPTIDE_QVALUE_COL:
    if ((scores_computed[LOGP_QVALUE_WEIBULL_XCORR] == true) &&
        (match->best_per_peptide == true)) {
      double qvalue = getScore(LOGP_PEPTIDE_QVALUE_WEIBULL);
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, qvalue);
    }
    break;
#endif
  case QVALUE_TDC_COL:
//    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(QVALUE_TDC));
//    }
    break;
  case QVALUE_MIXMAX_COL:
//    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(QVALUE_MIXMAX));
//    }
    break;
  case DECOY_XCORR_QVALUE_COL:
    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(DECOY_XCORR_QVALUE));
    }
    break;
  case DECOY_XCORR_PEP_COL:
    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       getScore(DECOY_XCORR_PEP));
    }
    break;
  case DECOY_EVALUE_QVALUE_COL:
    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(DECOY_EVALUE_QVALUE));
    }
    break;
  case DECOY_EVALUE_PEP_COL:
    if (null_peptide_ == false) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       getScore(DECOY_EVALUE_PEP));
    }
    break;
  
#ifdef NEW_COLUMNS
  case DECOY_XCORR_PEPTIDE_QVALUE_COL:
    if ( (null_peptide_ == false) && (best_per_peptide_ == true)) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
              getScore(DECOY_XCORR_PEPTIDE_QVALUE));
    }
    break;
#endif
  case PERCOLATOR_SCORE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(PERCOLATOR_SCORE));
    break;
  case PERCOLATOR_RANK_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getRank(PERCOLATOR_SCORE));
    break;
  case PERCOLATOR_QVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(PERCOLATOR_QVALUE));
    break;
  case PERCOLATOR_PEP_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(PERCOLATOR_PEP));
    break;
#ifdef NEW_COLUMNS
  case PERCOLATOR_PEPTIDE_QVALUE_COL:
    if ( match->best_per_peptide == true) {
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(PERCOLATOR_PEPTIDE_QVALUE));
    }
    break;
#endif
  case QRANKER_SCORE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(QRANKER_SCORE));
    break;
  case QRANKER_QVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(QRANKER_QVALUE));
    break;
  case QRANKER_PEP_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(QRANKER_PEP));
    break;
#ifdef NEW_COLUMNS
  case QRANKER_PEPTIDE_QVALUE_COL:
    if (match->best_per_peptide == true) {
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
              getScore(QRANKER_PEPTIDE_QVALUE));
    }
    break;
#endif
  case BY_IONS_MATCHED_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, b_y_matched);
    break;
  case BY_IONS_TOTAL_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, b_y_total);
    break;
  case DISTINCT_MATCHES_SPECTRUM_COL:
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       num_target_matches);
    break;
  case MATCHES_SPECTRUM_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
      num_target_matches);
    break;
  case DECOY_MATCHES_SPECTRUM_COL:
    if( null_peptide_ ){
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       num_decoy_matches);
    }
    break;
  case SEQUENCE_COL:
    {
      // this should get the sequence from the match, not the peptide
      char* sequence = getModSequenceStrWithMasses(
                get_mass_format_type_parameter("mod-mass-format"));
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, sequence);
      free(sequence);
    }
    break;
  case CLEAVAGE_TYPE_COL:
    {
      ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
      char* enzyme_string = enzyme_type_to_string(enzyme);
      DIGEST_T digestion = get_digest_type_parameter("digestion");
      char* digestion_string = digest_type_to_string(digestion);
      string cleavage_str = enzyme_string;
      cleavage_str += "-";
      cleavage_str += digestion_string;
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       cleavage_str.c_str() );
      free(enzyme_string);
      free(digestion_string);
    }
    break;
  case PROTEIN_ID_COL:
    {
      Peptide* peptide = getPeptide();
      string protein_ids_string = peptide->getProteinIdsLocations();
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       protein_ids_string.c_str());
    }
    break;
  case FLANKING_AA_COL:
    {
      Peptide* peptide = getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       flanking_aas);
      free(flanking_aas);
    }
    break;
  case ORIGINAL_TARGET_SEQUENCE_COL:
    if (null_peptide_ == true || OutputFiles::isConcat()) {
      char* seq = peptide_->getUnshuffledSequence();
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, seq);
      free(seq);
    }
    break;
  case ETA_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     collection->getCalibrationEta());
    break;
  case BETA_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     collection->getCalibrationBeta());
    break;
  case SHIFT_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     collection->getCalibrationShift());
    break;
  case CORR_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     collection->getCalibrationCorr());
    break;
  case INDEX_NAME_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
      database_index_name_);
    break;
    // values only for spectral-counts
  case RAW_SCORE_COL:
  case SIN_SCORE_COL:
  case NSAF_SCORE_COL:
  case DNSAF_SCORE_COL:
  case EMPAI_SCORE_COL:
  case PARSIMONY_RANK_COL:
    return;
  case NUMBER_MATCH_COLUMNS:
  case INVALID_COL:
    carp(CARP_FATAL, "Error in printing code (match.cpp).");
    break;
  }
}

/**
 * \brief Print the match information in tab delimited format to the given file
 *
 */
void Match::printTab(
  MatchCollection* collection,  ///< collection holding this match -in 
  MatchFileWriter* output_file,            ///< output stream -out
  Spectrum* spectrum,
  int num_target_matches,     ///< num matches in spectrum -in
  int num_decoy_matches ///< target matches for same spectrum -in
  ){

  carp(CARP_DETAILED_DEBUG, "Match::printTab: begin.");

  // Usually because no decoy file to print to.
  if( output_file == NULL ){ 
    return;
  }

  int b_y_total = getBYIonPossible();
  int b_y_matched = getBYIonMatched();
  
  // Print tab delimited fields
  int column_idx;
  for (column_idx = 0; column_idx < NUMBER_MATCH_COLUMNS; column_idx++) {
    carp(CARP_DETAILED_DEBUG,"print col:%i",column_idx);
    carp(CARP_DETAILED_DEBUG, "%s", get_column_header(column_idx));
    printOneMatchField(column_idx, 
                       collection,
                       output_file,
                       spectrum,
                       num_target_matches,
                       num_decoy_matches,
                       b_y_total,
                       b_y_matched);
  }
  output_file->writeRow();
  carp(CARP_DETAILED_DEBUG, "Match::printTab done.");
}

/**
 * shuffle the matches in the array between index start and end-1
 */
void Match::shuffleMatches(
  Match** match_array, ///< the match array to shuffle  
  int start_index,       ///< index of first element to shuffle
  int end_index          ///< index AFTER the last element to shuffle
  ){
  if( match_array == NULL ){
    carp(CARP_ERROR, "Cannot shuffle null match array.");
    return;
  }
  //  srandom(time(NULL));

  int match_idx = 0;
  for(match_idx = start_index; match_idx < end_index-1; match_idx++){
    Match* cur_match = match_array[match_idx];

    // pick a random index between match_index and end_index-1
    int rand_idx = get_random_number_interval(match_idx+1, end_index-1);

    //    fprintf(stderr, "%i values between %i and %i, rand %.4f, index %i\n",
    //            range, match_idx, end_index, rand_scaler, rand_idx);
    match_array[match_idx] = match_array[rand_idx];    
    match_array[rand_idx] = cur_match;
  }
}

/**
 *
 *\returns a match object that is parsed from the tab-delimited result file
 */
Match* Match::parseTabDelimited(
  MatchFileReader& result_file,  ///< the result file to parse PSMs -in
  Database* database, ///< the database to which the peptides are created -in
  Database* decoy_database ///< database with decoy peptides
  ) {

  string decoy_prefix = get_string_parameter("decoy-prefix");

  Match* match = new Match();

  Spectrum* spectrum = NULL;

  // this is a post_process match object
  match->post_process_match_ = true;

  Peptide* peptide = Peptide::parseTabDelimited(result_file, 
                                                database, 
                                                decoy_database);
  string index_name = result_file.getString(INDEX_NAME_COL);
  match->setDatabaseIndexName(index_name);

  if(peptide == NULL){
    carp(CARP_ERROR, "Failed to parse peptide (tab delimited)");
    // FIXME should this exit or return null. I think sometimes we can get
    // no peptides, which is valid, in which case NULL makes sense.
    // maybe this should be fixed at the output match level however.
    return NULL;
  }

  if ((result_file.empty(SP_SCORE_COL)) || (result_file.empty(SP_RANK_COL))){
    match -> match_scores_[SP] = NOT_SCORED;
    match -> match_rank_[SP] = 0;
  } else {
    match -> match_scores_[SP] = result_file.getFloat(SP_SCORE_COL);
    match -> match_rank_[SP] = result_file.getInteger(SP_RANK_COL);
  }

  match -> match_scores_[XCORR] = result_file.getFloat(XCORR_SCORE_COL);
  match -> match_rank_[XCORR] = result_file.getInteger(XCORR_RANK_COL);

  if (!result_file.empty(DELTA_CN_COL)) {
    match -> match_scores_[DELTA_CN] = result_file.getFloat(DELTA_CN_COL);
  }
  if (!result_file.empty(DELTA_LCN_COL)) {
    match -> match_scores_[DELTA_LCN] = result_file.getFloat(DELTA_LCN_COL);
  }

  if (!result_file.empty(EXACT_PVALUE_COL)){
    match -> match_scores_[TIDE_SEARCH_EXACT_PVAL] = result_file.getFloat(EXACT_PVALUE_COL);
    match -> match_scores_[TIDE_SEARCH_REFACTORED_XCORR] = result_file.getFloat(REFACTORED_SCORE_COL);
  }
  if (!result_file.empty(DECOY_XCORR_QVALUE_COL)){
    match->match_scores_[DECOY_XCORR_QVALUE] = result_file.getFloat(DECOY_XCORR_QVALUE_COL);
  }
  /* TODO I personally would like access to the raw p-value as well as the bonferonni corrected one (SJM).
  match -> match_scores[LOGP_WEIBULL_XCORR] = result_file.getFloat("logp weibull xcorr");
  */
  if (!result_file.empty(PVALUE_COL)){
    match->match_scores_[LOGP_BONF_WEIBULL_XCORR] = -log(result_file.getFloat(PVALUE_COL));
  }
  if (!result_file.empty(EVALUE_COL)){
    match->match_scores_[EVALUE]  = result_file.getFloat(EVALUE_COL);
  }
  if (!result_file.empty(PERCOLATOR_QVALUE_COL)){
    match->match_scores_[PERCOLATOR_QVALUE] = result_file.getFloat(PERCOLATOR_QVALUE_COL);
  }
  if (!result_file.empty(PERCOLATOR_SCORE_COL)){
    match->match_scores_[PERCOLATOR_SCORE] = result_file.getFloat(PERCOLATOR_SCORE_COL);
    match->match_rank_[PERCOLATOR_SCORE] = result_file.getInteger(PERCOLATOR_RANK_COL);
  }
  if (!result_file.empty(WEIBULL_QVALUE_COL)){
    match->match_scores_[LOGP_QVALUE_WEIBULL_XCORR] = result_file.getFloat(WEIBULL_QVALUE_COL);
  }
  if (!result_file.empty(QRANKER_SCORE_COL)){
    match->match_scores_[QRANKER_SCORE] = result_file.getFloat(QRANKER_SCORE_COL);
    match->match_scores_[QRANKER_QVALUE] = result_file.getFloat(QRANKER_QVALUE_COL);
  }

  if (!result_file.empty(BARISTA_SCORE_COL)) {
    match->match_scores_[BARISTA_SCORE] = result_file.getFloat(BARISTA_SCORE_COL);
    match->match_scores_[BARISTA_QVALUE] = result_file.getFloat(BARISTA_QVALUE_COL);
  }

  // get experiment size
  match->num_target_matches_ = 0;
  if (!result_file.empty(DISTINCT_MATCHES_SPECTRUM_COL)) {

    match->num_target_matches_ = result_file.getInteger(DISTINCT_MATCHES_SPECTRUM_COL);
  } else if (!result_file.empty(MATCHES_SPECTRUM_COL)) {
    match->num_target_matches_ = result_file.getInteger(MATCHES_SPECTRUM_COL);
  }
  if (match->num_target_matches_ == 0) {
    carp_once(CARP_WARNING, "num target matches=0, suppressing warning");
    match->ln_experiment_size_ = 0;
  } else {
    match->ln_experiment_size_ = log((FLOAT_T) match->num_target_matches_);
  }
  if (!result_file.empty(DECOY_MATCHES_SPECTRUM_COL)){
        match->num_decoy_matches_ = result_file.getInteger(DECOY_MATCHES_SPECTRUM_COL);
  }

  // parse spectrum
  if((spectrum = Spectrum::parseTabDelimited(result_file))== NULL){
    carp(CARP_ERROR, "Failed to parse spectrum (tab delimited).");
  }

  // spectrum specific features
  if (result_file.empty(BY_IONS_MATCHED_COL)){ 
    match -> b_y_ion_matched_ = 0;
  } else {
    match -> b_y_ion_matched_ = result_file.getInteger(BY_IONS_MATCHED_COL);
  }
  if (result_file.empty(BY_IONS_TOTAL_COL)) {
      match -> b_y_ion_possible_ = 0;
      match -> b_y_ion_fraction_matched_ = 0.0;
  } else {
    match -> b_y_ion_possible_ = result_file.getInteger(BY_IONS_TOTAL_COL);
    match -> b_y_ion_fraction_matched_ = 
      (FLOAT_T)match -> b_y_ion_matched_ /
      (FLOAT_T)match -> b_y_ion_possible_;
  }
  //parse match overall digestion
  match -> digest_ = string_to_digest_type((char*)result_file.getString(CLEAVAGE_TYPE_COL).c_str()); 

  if (!result_file.empty(PROTEIN_ID_COL) && 
    result_file.getString(PROTEIN_ID_COL).find(decoy_prefix) != string::npos) {
    match -> null_peptide_ = true;
  }

  //assign fields
  match -> peptide_sequence_ = NULL;
  match -> spectrum_ = spectrum;
  match -> peptide_ = peptide;
  
  return match;
}



/****************************
 * match get, set methods
 ***************************/


/**
 * Returns a heap allocated peptide sequence of the PSM
 * User must free the sequence.
 *
 * Go to top README for N,C terminus tryptic feature info.
 *
 *\returns the match peptide sequence, returns NULL if no sequence avaliable
 */
char* Match::getSequence() {
  // if post_process_match and has a null peptide you can't get sequence
  if(post_process_match_ && null_peptide_){
    carp(CARP_ERROR, 
         "Cannot retrieve null peptide sequence for post_process_match");
    return NULL;
  }
  
  if(peptide_sequence_ == NULL){
    peptide_sequence_ = peptide_->getSequence();
  }
  return my_copy_string(peptide_sequence_); 
}

/**
 * Returns a heap allocated peptide sequence of the PSM formatted with
 * the flanking amino acids and modifiation symbols.
 *
 * Sequence is in the form of X.SEQ.X where X is the flanking amino
 * acid or - if peptide is at the end of the protein.
 * Sequence may not be the same as for the peptide if this is for a
 * decoy database.
 *\returns The sqt-formatted peptide sequence for this match.
 */
char* Match::getSequenceSqt(){

  // get_match_mod_sequence (use method in case match->mod_seq == NULL) 
  MODIFIED_AA_T* mod_seq = getModSequence();
  if( mod_seq == NULL ){
    return NULL;
  }
  int length = getPeptide()->getLength();

  // turn it into string
  char* seq = modified_aa_string_to_string_with_masses(mod_seq, length,
                get_mass_format_type_parameter("mod-mass-format"));

  // get peptide flanking residues 
  char c_term = peptide_->getCTermFlankingAA();
  char n_term = peptide_->getNTermFlankingAA();

  // allocate seq + 4 length array
  char* final_string = (char*)mycalloc((strlen(seq)+5), sizeof(char));

  // copy pieces in
  final_string[0] = n_term;
  final_string[1] = '.';
  strcpy(&final_string[2], seq);
  final_string[strlen(seq) + 2] = '.';
  final_string[strlen(seq) + 3] = c_term;
  final_string[strlen(seq) + 4] = '\0';

  carp(CARP_DETAILED_DEBUG, "start string %s, final %s", seq, final_string);

  // delete mod seq and string version
  free(seq);
  free(mod_seq);
  return final_string;
}

/**
 * \brief Returns a newly allocated modified_aa sequence of the PSM.
 * Sequence is the same as the peptide, if target match or is a
 * shuffled sequence if a null (decoy) match.  If match field
 * 'mod_sequence' is non NULL, returns a copy of that value, otherwise
 * fills that field and returns a copy of the value.
 *
 * \returns the match peptide sequence, returns NULL if no sequence avaliable
 */
MODIFIED_AA_T* Match::getModSequence()
{

  // if post_process_match and has a null peptide you can't get sequence
  if(post_process_match_ && null_peptide_){
    return NULL;
  }

  if(mod_sequence_ == NULL){
    mod_sequence_ = peptide_->getModifiedAASequence();
  }

  return copy_mod_aa_seq(mod_sequence_, peptide_->getLength());
}

/**
 * \brief Returns a newly allocated string of sequence including any
 * modifications represented as symbols (*,@,#, etc) following the
 * modified residue. 
 * \returns The peptide sequence of the match including modification
 * characters. 
 */
char* Match::getModSequenceStrWithSymbols(){

  // if post_process_match and has a null peptide you can't get sequence
  if(post_process_match_ && null_peptide_){
    return NULL;
  }

  if(mod_sequence_ == NULL){
    mod_sequence_ = peptide_->getModifiedAASequence();
  }

  return modified_aa_string_to_string_with_symbols(mod_sequence_, 
                                      peptide_->getLength());
}

/**
 * \brief Returns a newly allocated string of sequence including any
 * modifications represented as mass values in brackets following the
 * modified residue. If merge_masses is true, the sum of multiple
 * modifications on one residue are printed.  If false, each mass is
 * printed in a comma-separated list.
 * \returns The peptide sequence of the match including modification
 * masses. 
 */
char* Match::getModSequenceStrWithMasses(
 MASS_FORMAT_T mass_format)
{

  // if post_process_match and has a null peptide you can't get sequence
  if(post_process_match_ && null_peptide_){
    return NULL;
  }

  if(mod_sequence_ == NULL){
    mod_sequence_ = peptide_->getModifiedAASequence();
  }

  return modified_aa_string_to_string_with_masses(mod_sequence_, 
                                                  peptide_->getLength(),
                                                  mass_format);
}
/**
 * Must ask for score that has been computed
 *\returns the match_mode score in the match object
 */
FLOAT_T Match::getScore (
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  ) const
{
  return match_scores_[match_mode];
}

/**
 * sets the match score
 */
void Match::setScore(
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  FLOAT_T match_score ///< the score of the match -in
  )
{
  match_scores_[match_mode] = match_score;
}


/**
 * set the custom match score
 */
void Match::setCustomScore(
  const std::string& match_score_name, ///< the name of the score -in
  FLOAT_T match_score ///< the score of the match -in
  ) {

  match_custom_scores_[match_score_name] = match_score;

}

/**
 * get the custom score
 */
bool Match::getCustomScore(
  const std::string& match_score_name, ///< the name of the score -in
  FLOAT_T& score ///< the value of the score -out
  ) {

  if (match_custom_scores_.find(match_score_name) == match_custom_scores_.end()) {
    carp(CARP_ERROR, "custom match score:%s doesn't exist!", match_score_name.c_str());
    return false;
  }

  score = match_custom_scores_[match_score_name];
  return true;

}

void Match::getCustomScoreNames(
  vector<string>& custom_score_names
  ) {

  custom_score_names.clear();

  for (map<string,FLOAT_T>::iterator iter = match_custom_scores_.begin();
    iter != match_custom_scores_.end();
      ++iter) {
    custom_score_names.push_back(iter->first);
  }

}

/**
 * sets the file index for this match
 */
void Match::setFileIndex(
  int file_idx ///< file index to set
  ) {

  file_idx_ = file_idx;
}

/**
 * \returns the file index for this match
 */
int Match::getFileIndex() {
  return(file_idx_);
}

int Match::findFileIndex(const string& file_path, bool match_stem) {
  for (size_t idx = 0; idx < file_paths_.size(); idx++) {
    if (file_path == file_paths_[idx]) {
      return idx;
    }
  }
  if (match_stem) {
    for (size_t idx = 0; idx < file_paths_.size(); idx++) {
      boost::filesystem::path boost_path(file_paths_[idx]);
      if (file_path == boost_path.stem().string()) {
        return idx;
      }
    }
  }
  return -1;
}

/**
 * sets the file path for this match
 * \returns the associated file index
 */
int Match::setFilePath(
  const string& file_path ///< file path to set
  ) {
  file_idx_ = findFileIndex(file_path);
  if (file_idx_ == -1) {
    file_idx_ = file_paths_.size();
    file_paths_.push_back(file_path);
  }
  return file_idx_;
}

/**
 * \returns the file path for this match
 */ 
string Match::getFilePath() {
  return(getFilePath(file_idx_));
}

string Match::getFilePath(int file_idx) {
  if (file_idx == -1 || file_paths_.empty()) {
    return string("");
  } else {
    return(file_paths_[file_idx]);
  }
}



bool Match::isDecoy() {

  return null_peptide_;

}

/**
 * sets whether the match is post process or not
 */
void Match::setPostProcess(
  bool post_process ///< whether the match is post process or not
) {
  post_process_match_ = post_process;
}

/**
 * Must ask for score that has been computed
 *\returns the match_mode rank in the match object
 */
int Match::getRank(
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match_rank_[match_mode];
}

/**
 * sets the rank of the match
 */
void Match::setRank(
  SCORER_TYPE_T match_mode, ///< the working mode -in
  int match_rank ///< the rank of the match -in
  )
{
  match_rank_[match_mode] = match_rank;
}

/**
 *\returns the spectrum in the match object
 */
Spectrum* Match::getSpectrum()
{
  return spectrum_;
}

/**
 *\returns the peptide in the match object
 */
Peptide* Match::getPeptide()
{
  return peptide_;
}

/**
 * sets whether the match is a null peptide match or not
 */
void Match::setNullPeptide(
  bool null_peptide ///< whether the match is a null peptide match or not
) {
  null_peptide_ = null_peptide;
}

/**
 * gets the match if it is a null_peptide match
 *\returns true if match is null peptide, else false
 */
bool Match::getNullPeptide()
{
  return null_peptide_;
}

void Match::setZState(
  SpectrumZState& zstate) {

  zstate_ = zstate;

}

SpectrumZState& Match::getZState() {

  return zstate_;

}

/**
 * gets the match charge
 */
int Match::getCharge()
{
  return getZState().getCharge();
}

/**
 * /returns the match neutral mass
 */
FLOAT_T Match::getNeutralMass()
{
  return getZState().getNeutralMass();
}

/**
 * sets the match ln_experiment_size
 */
void Match::setLnExperimentSize(
  FLOAT_T ln_experiment_size ///< the ln_experiment_size value of PSM -in
  )
{
  ln_experiment_size_ = ln_experiment_size;
}

/**
 * gets the match ln_experiment_size
 */
FLOAT_T Match::getLnExperimentSize()
{
  return ln_experiment_size_;
}

/**
 * Sets the total number of target matches searched for this spectrum.
 */
void Match::setTargetExperimentSize(int num_matches){
  num_target_matches_ = num_matches;
}

/**
 * \returns The total number of target matches searched for this spectrum.
 */
int Match::getTargetExperimentSize(){
  return num_target_matches_;
}

/**
 * \returns The total number of decoy matches searched for this
 * spectrum if this is a match to a decoy spectrum.
 */
int Match::getDecoyExperimentSize(){
  return num_decoy_matches_;
}

/**
 * sets the match b_y_ion information
 */
void Match::setBYIonInfo(
  Scorer* scorer ///< the scorer from which to extract information -in
  ) {

  b_y_ion_fraction_matched_ = scorer->getSpBYIonFractionMatched(); 
  b_y_ion_matched_ = scorer->getSpBYIonMatched(); 
  b_y_ion_possible_ = scorer->getSpBYIonPossible(); 
}

void Match::setBYIonFractionMatched(
  FLOAT_T fraction_matched
  ) {

  b_y_ion_fraction_matched_ = fraction_matched;
}

void Match::calcBYIonFractionMatched() {
  b_y_ion_fraction_matched_ = (FLOAT_T)b_y_ion_matched_ / (FLOAT_T)b_y_ion_possible_;
}

/**
 * gets the match b_y_ion_fraction_matched
 */
FLOAT_T Match::getBYIonFractionMatched()
{
  return b_y_ion_fraction_matched_;
}

/**
 * sets the match b_y_ion_matched
 */
void Match::setBYIonMatched(int matched) {
  b_y_ion_matched_ = matched;
  if (b_y_ion_matched_ > 0 && b_y_ion_possible_ > 0) {
    calcBYIonFractionMatched();
  }
  
}

/**
 * gets the match b_y_ion_matched
 */
int Match::getBYIonMatched()
{
  return b_y_ion_matched_;
}


/**
 * sets the match b_y_ion_possible
 */
void Match::setBYIonPossible(int possible) {
  b_y_ion_possible_ = possible;
  if (b_y_ion_possible_ > 0 && b_y_ion_matched_ > 0) {
    calcBYIonFractionMatched();
  }
} 

/**
 * gets the match b_y_ion_possible
 */
int Match::getBYIonPossible()
{
  return b_y_ion_possible_;
}

/**
 *Increments the pointer count to the match object
 */
void Match::incrementPointerCount()
{
  ++pointer_count_;
}


/**
 * Set the best-per-peptide Boolean to true.
 */
void Match::setBestPerPeptide() {
  best_per_peptide_ = true;
}


/************************************************
 * TODO: Why are these here?
 ************************************************/

/**
 * \brief prints both variable and static modifications for 
 *  peptide sequence in xml format to the specificed output file
 *
 *
 */
void print_modifications_xml(
  const char* mod_seq,
  const char* pep_seq,
  FILE* output_file
){
  map<int, double> var_mods;
  map<int, double> static_mods;
  carp(CARP_DEBUG,"print_modifications_xml:%s %s", mod_seq, pep_seq);
  // variable modifications
  int mod_precision = get_int_parameter("mod-precision");
  find_variable_modifications(var_mods, mod_seq);
  if (!var_mods.empty()){
    fprintf(output_file, 
            "<modification_info modified_peptide=\"%s\">\n",
            mod_seq);
   carp(CARP_DEBUG,
            "<modification_info modified_peptide=\"%s\">\n",
            mod_seq);
        for (map<int, double>::iterator it = var_mods.begin()
           ; it != var_mods.end(); ++it){
      fprintf(output_file, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              (*it).first,   //index
              mod_precision, (*it).second); //mass
      carp(CARP_DEBUG, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              (*it).first,   //index                                                                                                                                                                                                            
              mod_precision, (*it).second); //mass           
    }
    fprintf(output_file, "</modification_info>\n");
  }

  // static modifications
  find_static_modifications(static_mods, var_mods, pep_seq);
  if (!static_mods.empty()){
    carp(CARP_DEBUG, "<modification_info modified_peptide=\"%s\">\n",
            pep_seq);
    fprintf(output_file, "<modification_info modified_peptide=\"%s\">\n",
            pep_seq);
    for (map<int, double>::iterator it = static_mods.begin(); 
         it != static_mods.end(); ++it){
      fprintf(output_file, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              (*it).first,   //index
              mod_precision, (*it).second); //mass
      carp(CARP_DEBUG, "<mod_aminoacid_mass position=\"%i\" mass=\"%.*f\"/>\n",
              (*it).first,   //index                                                                                                                                                                                                            
              mod_precision, (*it).second); //mass 
    }
    fprintf(output_file, "</modification_info>\n");
  }

  static_mods.clear();
  var_mods.clear();

}


/**
 * \brief takes an empty mapping of index to mass
 * and extract information from mod sequence fill
 * up map
 */
void find_variable_modifications(
 map<int, double>& mods,
 const char* mod_seq
){
  
    int seq_index = 1;
    const char* amino = mod_seq;
    const char* end = NULL;
    const char* start = NULL;
    // Parse returned string to find modifications within
    // brackets
    while (*(amino) != '\0' && *(amino+1) != '\0'){
      if (*(amino+1) =='['){
        start = amino+2;
        end = amino+2;
        while (*end != ']'){
          end++;
        }
        char* mass  = (char *) mymalloc(sizeof(char)*(end-start+1));
        strncpy(mass, start, end-start);
        mass[end-start] = '\0';
        mods[seq_index] = atof(mass);
        amino = end;
        free(mass);
      } else if (*(amino+1) < 'A' || *(amino+1) > 'Z'){ // a mod symbol
        double mass = 0; // sum up all adjacent symbols
        end = amino + 1;
        while( *end < 'A' || *end > 'Z' ){
          mass += get_mod_mass_from_symbol(*end);
          end++;
        }
        mods[seq_index] = mass;
      }
      seq_index++;
      amino++;
    }
}

/**
 * \brief takes an empty mapping of index to mass
 * of static mods and a full mapping of var mods
 * to fill up the mapping of static mods
 */
void find_static_modifications(
  map<int, double>& static_mods,
  map<int, double>& var_mods,
  const char* peptide_sequence
){
  const char* seq_iter = peptide_sequence;
  
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");
  
  char aa[2];
  aa[1] = '\0';
  int seq_index = 1;
  while ((*seq_iter) != '\0'){
    *aa = (*seq_iter);
    // write static mod if user requested static mod and also
    // if there is no variable mod on the same character
    if (get_double_parameter( (const char *)aa)!= 0 && 
        var_mods.find(seq_index) == var_mods.end()){

      //double mass = get_mass_amino_acid(*seq_iter, isotopic_type);
      static_mods[seq_index] = get_double_parameter( (const char*)aa);
    }
    seq_iter++;
    seq_index++;
  }

}


/**
 * \brief Counts the number of internal cleavages
 *
 */
int get_num_internal_cleavage(const char* peptide_sequence, ENZYME_T enzyme){
  // get number of internal cleavages
  int num_missed_cleavages = 0;
  const char* seq_iter = peptide_sequence;
  
  while (*(seq_iter+1) != '\0'){
    if (ProteinPeptideIterator::validCleavagePosition(seq_iter, enzyme)){
      num_missed_cleavages++;
    }
    seq_iter++;
  }
  return num_missed_cleavages;
}

/**
 * \brief Returns whether the nterm and cterm of a peptide are proper cleavages
 */
void get_terminal_cleavages(
  const char* peptide_sequence, ///< peptide sequenc
  const char flanking_aas_prev, ///< amino acid before cleavage (n-term)
  const char flanking_aas_next, ///< amino acid after cleavage (c-term)
  ENZYME_T enzyme, ///< Enzyme used in cleavage
  bool& nterm, ///< -out is nterminus from a proper cleavage
  bool& cterm ///< -out is cterminus from a proper cleavage?
  ) {
  int num_tol_term = 0;
  char cleavage[3];
  cleavage[2] = '\0';
  cleavage[0] = flanking_aas_prev;
  cleavage[1] = peptide_sequence[0];
  nterm = flanking_aas_prev == '-' || ProteinPeptideIterator::validCleavagePosition(cleavage, enzyme);

  cleavage[0] = peptide_sequence[strlen(peptide_sequence)-1];
  cleavage[1] = flanking_aas_next;
  cterm = flanking_aas_next == '-' || ProteinPeptideIterator::validCleavagePosition(cleavage, enzyme);
}

/**
 * \brief Counts the number of terminal cleavage. Either 0, 1, or 2
 *
 */
int get_num_terminal_cleavage(
  const char* peptide_sequence, 
  const char flanking_aas_prev,
  const char flanking_aas_next,
  ENZYME_T enzyme
  ){

  bool cterm, nterm;

  get_terminal_cleavages(peptide_sequence, 
    flanking_aas_prev, 
    flanking_aas_next, 
    enzyme,
    nterm, cterm);

  int num_tol_term = 0;
  if (nterm) { 
    num_tol_term++;
  }
  if (cterm) { 
    num_tol_term++;
  }
  return num_tol_term;
}

void Match::setDatabaseIndexName(string index_name){
  database_index_name_ = index_name;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

