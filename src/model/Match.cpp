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
#include "parameter.h"
#include "util/GlobalParams.h"
#include "Scorer.h" 
#include "Match.h" 
#include "MatchCollection.h" 
#include "Peptide.h"
#include "util/Params.h"
#include "util/FileUtils.h"

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
  // default is not a null peptide match
  null_peptide_ = false;
  peptide_sequence_ = NULL;
  mod_sequence_ = NULL;
  post_process_match_ = 0;
  ln_experiment_size_ = 0;
  num_target_matches_ = 0;
  best_per_peptide_ = false;
  file_idx_ = -1;
  decoy_idx_ = -1;
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
             const SpectrumZState& zstate, ///< the charge/mass of the spectrum
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
  /*
  if (peptide_sequence_ != NULL){
    free(peptide_sequence_);
  }
  */
  if (mod_sequence_ != NULL){
    free(mod_sequence_);
  }
}

bool Match::ScoreComparer::operator() (const Match* x, const Match* y) {
  FLOAT_T scoreX = x->getScore(type_);
  FLOAT_T scoreY = y->getScore(type_);
  return less_ ? Match::ScoreLess(scoreX, scoreY) : Match::ScoreGreater(scoreX, scoreY);
}

bool Match::ScoreLess(FLOAT_T x, FLOAT_T y) {
  if (isnan(x) || isinf(x)) {
    x = std::numeric_limits<FLOAT_T>::max();
  }
  if (isnan(y) || isinf(y)) {
    y = std::numeric_limits<FLOAT_T>::max();
  }
  return x < y;
}

bool Match::ScoreGreater(FLOAT_T x, FLOAT_T y) {
  if (isnan(x) || isinf(x)) {
    x = -std::numeric_limits<FLOAT_T>::max();
  }
  if (isnan(y) || isinf(y)) {
    y = -std::numeric_limits<FLOAT_T>::max();
  }
  return x > y;
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
  if (!file) {
    carp(CARP_ERROR, "Cannot print match to sqt file from null input");
    return;
  }

  Peptide* peptide = getPeptide();
  // this should get the sequence from the match, not the peptide
  char* sequence = getSequenceSqt();

  int b_y_total = getScore(BY_IONS_TOTAL);
  if (b_y_total == NOT_SCORED) {
    b_y_total = 0;
  }
  int b_y_matched = getScore(BY_IONS_MATCHED);
  if (b_y_matched == NOT_SCORED) {
    b_y_matched = 0;
  }
  
  FLOAT_T delta_cn = getScore(DELTA_CN);
  FLOAT_T score_main = getScore(XCORR);

  // write format string with variable precision
  int precision = Params::GetInt("precision");

  // print match info
  if (exact_pval_search_) {
    fprintf(file, "M\t%i\t%i\t%.*f\t%.2f\t%.*g\t%.*g\t%.*g\t%i\t%i\t%s\tU\n",
            getRank(XCORR),
            getRank(SP),
            Params::GetInt("mass-precision"),
            peptide->calcModifiedMass() + MASS_PROTON,
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
            Params::GetInt("mass-precision"),
            peptide->calcModifiedMass() + MASS_PROTON,
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
  string protein_id;
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
      rand = Params::GetString("decoy-prefix"); 
    }

    // print match info (locus line), add "decoy-prefix" to locus name for decoys
    fprintf(file, "L\t%s%s\n", rand.c_str(), protein_id.c_str());      
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
      double peptide_mass = peptide->calcModifiedMass();
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
  case RESIDUE_EVIDENCE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getScore(RESIDUE_EVIDENCE_SCORE));
    break;
  case BOTH_PVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getScore(BOTH_PVALUE));
    break;
  case SIDAK_ADJUSTED_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                     getScore(SIDAK_ADJUSTED));
    break;
  case XCORR_RANK_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getRank(!Params::GetBool("exact-p-value") ? XCORR : TIDE_SEARCH_EXACT_PVAL));
    break;
  case RESIDUE_RANK_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getRank(RESIDUE_EVIDENCE_SCORE));
    break;
  case BOTH_PVALUE_RANK:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                     getRank(BOTH_PVALUE));
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
      } else {
        double p_value = log_pvalue != numeric_limits<FLOAT_T>::infinity() ?
          exp(-1 * log_pvalue) : 0;
        output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                         p_value);
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
  case MODIFICATIONS_COL:
    {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       getPeptide()->getModsString());
    }
    break;
  case CLEAVAGE_TYPE_COL:
    {
      ENZYME_T enzyme = GlobalParams::getEnzyme();
      const char* enzyme_string = enzyme_type_to_string(enzyme);
      DIGEST_T digestion = GlobalParams::getDigestion();
      const char* digestion_string = digest_type_to_string(digestion);
      string cleavage_str = enzyme_string;
      cleavage_str += "-";
      cleavage_str += digestion_string;
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
                                       cleavage_str.c_str() );
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
  case TARGET_DECOY_COL:
    if ( null_peptide_ ) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       "decoy");
    } else {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       "target");
    }      
    break;
  case ORIGINAL_TARGET_SEQUENCE_COL:
    if (null_peptide_ || Params::GetBool("concat")) {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                       peptide_->getUnshuffledSequence());
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

  int b_y_total = getScore(BY_IONS_TOTAL);
  if (b_y_total == NOT_SCORED) {
    b_y_total = 0;
  }
  int b_y_matched = getScore(BY_IONS_MATCHED);
  if (b_y_matched == NOT_SCORED) {
    b_y_matched = 0;
  }
  
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
      if (file_path == FileUtils::Stem(file_paths_[idx])) {
        return idx;
      }
    }
  }
  return -1;
}

int Match::addUniqueFilePath(const string& path, bool match_stem) {
  int idx = findFileIndex(path, match_stem);
  if (idx == -1) {
    idx = file_paths_.size();
    file_paths_.push_back(path);
    carp(CARP_INFO, "Assigning index %d to %s.", idx, path.c_str());
  }
  return idx;
}

/**
 * sets the file path for this match
 * \returns the associated file index
 */
int Match::setFilePath(
  const string& file_path ///< file path to set
  ) {
  file_idx_ = addUniqueFilePath(file_path);
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
    return "";
  } else {
    return(file_paths_[file_idx]);
  }
}

int Match::decoyIndex() const {
  return decoy_idx_;
}

void Match::setDecoyIndex(int value) {
  decoy_idx_ = value;
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
bool Match::getNullPeptide() {
  return null_peptide_;
}

void Match::setZState(SpectrumZState& zstate) {
  zstate_ = zstate;
}

SpectrumZState& Match::getZState() {
  return zstate_;
}

/**
 * gets the match charge
 */
int Match::getCharge() {
  return getZState().getCharge();
}

/**
 * /returns the match neutral mass
 */
FLOAT_T Match::getNeutralMass() {
  return getZState().getNeutralMass();
}

/**
 * sets the match ln_experiment_size
 */
void Match::setLnExperimentSize(
  FLOAT_T ln_experiment_size ///< the ln_experiment_size value of PSM -in
) {
  ln_experiment_size_ = ln_experiment_size;
}

/**
 * gets the match ln_experiment_size
 */
FLOAT_T Match::getLnExperimentSize() {
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
 *Increments the pointer count to the match object
 */
void Match::incrementPointerCount() {
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

