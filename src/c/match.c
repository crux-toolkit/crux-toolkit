/*************************************************************************//**
 * \file match.c
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * DESCRIPTION: Object for matching a peptide and a spectrum, generate
 * a preliminary score(e.g., Sp) 
 *
 * REVISION: $Revision: 1.78 $
 * REVISION: $Revision: 1.78 $
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "parse_arguments.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "ion.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "scorer.h" 
#include "match.h" 
#include "match_collection.h" 
#include "generate_peptides_iterator.h" 
#include "peptide.h"

/**
 * README!
 * <Issues on the overall_type field in match struct>
 * 
 * Outstanding question: How do you determine the
 * peptide trypticity for multiple protein sources?
 * 
 * For example, if one protein src is tryptic while the other is 
 * not tryptic, what is the peptide trypticity used for feature and 
 * shuffling peptide sequence?
 *
 * Currently, we use the "Tryptic wins all" approach, where
 * if N-terminus is tryptic in any of the src proteins we claim it
 * tryptic on the N terminus. Same applies for C-terminus.
 * Thus, even if no peptide is tryptic in any of its src protein,
 * if there are one src protein where it is N tryptic and another where
 * it is C tryptic, overall we will call the peptide in the match tryptic.
 * 
 * There are for sure other methods, for example to randomly sample the src 
 * protein and consider that protein as its src or to shuffle the flanking 
 * sequence of the peptide in each of peptide and randomly sample from the 
 * shuffled flanking sequence to determine the shuffled peptide's trypticity.
 *
 */

/**
 *\struct match
 *\brief An object that stores the score & rank for each pepide-spectrum match
 */
struct match{
  SPECTRUM_T* spectrum; ///< the spectrum we are scoring with
  PEPTIDE_T* peptide;  ///< the peptide we are scoring
  FLOAT_T match_scores[_SCORE_TYPE_NUM]; 
    ///< array of scores, one for each type (index with SCORER_TYPE_T) 
  int match_rank[_SCORE_TYPE_NUM];  
    ///< rank of this match for each type scored (index with SCORER_TYPE_T)
  int pointer_count; 
    ///< number of pointers to this match object (when reach 0, free memory)
  FLOAT_T b_y_ion_fraction_matched; 
    ///< the fraction of the b, y ion matched while scoring for SP
  int b_y_ion_matched; ///< number of b, y ion matched while scoring SP
  int b_y_ion_possible; ///< number of possible b, y ion while scoring SP
  BOOLEAN_T null_peptide; ///< Is the match a null (decoy) peptide match?
  char* peptide_sequence; ///< peptide sequence is that of peptide or shuffled
  MODIFIED_AA_T* mod_sequence; ///< seq of peptide or shuffled if null peptide
  DIGEST_T digest;
  //  PEPTIDE_TYPE_T overall_type; 
    ///< overall peptide trypticity, set in set_match_peptide, see README above
  int charge; ///< the charge state of the match 
  // post_process match object features
  // only valid when post_process_match is TRUE
  BOOLEAN_T post_process_match; ///< Is this a post process match object?
  FLOAT_T delta_cn; ///< the difference in top and second Xcorr scores
  FLOAT_T ln_delta_cn; ///< the natural log of delta_cn
  FLOAT_T ln_experiment_size; 
     ///< natural log of total number of candidate peptides evaluated
};

/**
 * \returns a new memory allocated match
 */
MATCH_T* new_match(void){
  MATCH_T* match = (MATCH_T*)mycalloc(1, sizeof(MATCH_T));
  
  // initialize score, rank !!!!DEBUG
  int index = 0;
  for(index = 0; index < _SCORE_TYPE_NUM; ++index){
    match->match_rank[index] = 0;
    match->match_scores[index] = NOT_SCORED;
  }
  
  ++match->pointer_count;

  // default is not a null peptide match
  match->null_peptide = FALSE;

  // set default as not tryptic
  // a full evaluation is done when set peptide
  //  match->overall_type = NOT_TRYPTIC;

  return match;
}

/**
 * free the memory allocated match
 * spectrum is not freed by match
 */
void free_match(
  MATCH_T* match ///< the match to free -in
  )
{
  --match->pointer_count;
  
  // only free match when pointer count reaches
  if(match->pointer_count == 0){

    // but aren't there multiple matches pointing to the same peptide?
    // if so, create a new free_shallow_match which doesn't touch the members
    if (match->peptide != NULL){
      free_peptide(match->peptide);
    }
    if(match->post_process_match && match->spectrum !=NULL){
      free_spectrum(match->spectrum);
    }
    if (match->peptide_sequence != NULL){
      free(match->peptide_sequence);
    }
    if (match->mod_sequence != NULL){
      free(match->mod_sequence);
    }

    free(match);  
  }
}

/* ************ SORTING FUNCTIONS ****************** */
/**
 * Compare two matches by spectrum scan number.
 * \return 0 if scan numbers are equal, -1 if a is smaller, 1 if b is
 * smaller.
 */
int compare_match_spectrum(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  ){

  SPECTRUM_T* spec_a = get_match_spectrum((*match_a));
  SPECTRUM_T* spec_b = get_match_spectrum((*match_b));
  int scan_a = get_spectrum_first_scan(spec_a);
  int scan_b = get_spectrum_first_scan(spec_b);
  int charge_a = get_match_charge((*match_a));
  int charge_b = get_match_charge((*match_b));

  if( scan_a < scan_b ){
    return -1;
  }else if( scan_a > scan_b ){
    return 1;
  }else{// else scan numbers equal
    if(charge_a < charge_b){
      return -1;
    }else if(charge_a > charge_b){
      return 1;
    }
  }// else scan number and charge equal
  return 0;
}

/**
 * compare two matches, used for qsort
 * \returns 0 if sp scores are equal, -1 if a is less than b, 1 if a
 * is greather than b.
 */
int compare_match_sp(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{
  // might have to worry about cases below 1 and -1
  // return(int)((*match_b)->match_scores[SP] - (*match_a)->match_scores[SP]);

  if((*match_b)->match_scores[SP] > (*match_a)->match_scores[SP]){
    return 1;
  }
  else if((*match_b)->match_scores[SP] < (*match_a)->match_scores[SP]){
    return -1;
  }
  return 0;

}

/**
 * Compare two matches by spectrum scan number and sp score, used for qsort.
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if sp score of match a is less than
 * match b.  1 if scan number and sp are equal, else 0.
 */
int compare_match_spectrum_sp(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  ){

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_sp(match_a, match_b);
  }

  return return_me;
}

/**
 * compare two matches, used for qsort
 * \returns 0 if xcorr scores are equal, -1 if a is less than b, 1 if a
 * is greather than b.
 */
int compare_match_xcorr(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{

  if((*match_b)->match_scores[XCORR] > (*match_a)->match_scores[XCORR]){
    return 1;
  }
  else if((*match_b)->match_scores[XCORR] < (*match_a)->match_scores[XCORR]){
    return -1;
  }
  return 0;

}

/**
 * Compare two matches by spectrum scan number and xcorr, used for qsort.
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_xcorr(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
){

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_xcorr(match_a, match_b);
  }

  return return_me;
}

/**
 * compare two matches, used for qsort
 * \returns the difference between p_value (LOGP_BONF_WEIBULL_XCORR)
 * score in match_a and match_b 
 */
int compare_match_p_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  ){

  if((*match_b)->match_scores[LOGP_BONF_WEIBULL_XCORR] 
     > (*match_a)->match_scores[LOGP_BONF_WEIBULL_XCORR]){
    return 1;
  }
  else if((*match_b)->match_scores[LOGP_BONF_WEIBULL_XCORR] 
          < (*match_a)->match_scores[LOGP_BONF_WEIBULL_XCORR]){
    return -1;
  }
  return 0;

}


/**
 * compare two matches, used for qsort
 * The smaller the Q value the better!!!, this is opposite to other scores
 * \returns 0 if q-value scores are equal, -1 if a is less than b, 1 if a
 * is greather than b.
 */
int compare_match_q_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{

  if((*match_b)->match_scores[Q_VALUE] < (*match_a)->match_scores[Q_VALUE]){
    return 1;
  }
  else if((*match_b)->match_scores[Q_VALUE] > (*match_a)->match_scores[Q_VALUE]){
    return -1;
  }
  return 0;
}

/**
 * compare two matches, used for qsort
 * The smaller the Q value the better!!!, this is opposite to other scores
 * \returns 0 if q-value scores are equal, -1 if a is less than b, 1 if a
 * is greather than b.
 */
int compare_match_qranker_q_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{

  if((*match_b)->match_scores[QRANKER_Q_VALUE] < 
      (*match_a)->match_scores[QRANKER_Q_VALUE]){
    return 1;
  }
  else if((*match_b)->match_scores[QRANKER_Q_VALUE] > 
          (*match_a)->match_scores[QRANKER_Q_VALUE]){
    return -1;
  }
  return 0;
}

/**
 * Compare two matches by spectrum scan number and q-value, used for qsort.
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_q_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
){

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_q_value(match_a, match_b);
  }

  return return_me;
}

/**
 * Compare two matches by spectrum scan number and qranker q-value, 
 * used for qsort.
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_qranker_q_value(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
){

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_q_value(match_a, match_b);
  }

  return return_me;
}


/**
 * compare two matches, used for PERCOLATOR_SCORE
 * \returns 0 if percolator scores are equal, -1 if a is less than b,
 * 1 if a is greather than b.
 */
int compare_match_percolator_score(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{
  if((*match_b)->match_scores[PERCOLATOR_SCORE] > (*match_a)->match_scores[PERCOLATOR_SCORE]){
    return 1;
  }
  else if((*match_b)->match_scores[PERCOLATOR_SCORE] < (*match_a)->match_scores[PERCOLATOR_SCORE]){
    return -1;
  }
  return 0;

}

/**
 * Compare two matches by spectrum scan number and percolator score,
 * used for qsort. 
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_percolator_score(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_percolator_score(match_a, match_b);
  }

  return return_me;
}

/**
 * compare two matches, used for QRANKER_SCORE
 * \returns 0 if qranker scores are equal, -1 if a is less than b,
 * 1 if a is greather than b.
 */
int compare_match_qranker_score(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
)
{
  if((*match_b)->match_scores[QRANKER_SCORE] > (*match_a)->match_scores[QRANKER_SCORE]){
    return 1;
  }
  else if((*match_b)->match_scores[QRANKER_SCORE] < (*match_a)->match_scores[QRANKER_SCORE]){
    return -1;
  }
  return 0;

}

/**
 * Compare two matches by spectrum scan number and qranker score,
 * used for qsort. 
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_qranker_score(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{

  int return_me = compare_match_spectrum( match_a, match_b );

  if( return_me == 0 ){
    return_me = compare_match_qranker_score(match_a, match_b);
  }

  return return_me;
}

/**
 * Compare two matches by spectrum scan number and q-value (from the decoys and xcorr score),
 * used for qsort. 
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_decoy_xcorr_qvalue(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{
  // delete this, just for the compiler
  if( match_a == NULL || match_b == NULL ){
    return 0;
  }
  carp(CARP_FATAL, "HEY, you haven't implemented sorting by decoy-qvalue yet!");
  return 0; // Return value to avoid compiler warning.
}

/**
 * Compare two matches by spectrum scan number and q-value (from the decoys and weibull est p-values),
 * used for qsort. 
 * \returns -1 if match a spectrum number is less than that of match b
 * or if scan number is same, if score of match a is less than
 * match b.  1 if scan number and score are equal, else 0.
 */
int compare_match_spectrum_decoy_pvalue_qvalue(
  MATCH_T** match_a, ///< the first match -in  
  MATCH_T** match_b  ///< the scond match -in
  )
{
  // delete this, just for the compiler
  if( match_a == NULL || match_b == NULL ){
    return 0;
  }
  carp(CARP_FATAL, "HEY, you haven't implemented sorting by decoy-qvalue yet!");
  return 0; // Return value to avoid compiler warning.
}

/* ****** End of sorting functions ************/

/**
 * FIXME this should be capable of returning with an error
 * print the information of the match
 */
void print_match(
  MATCH_T* match, ///< the match to print -in  
  FILE* file, ///< output stream -out
  BOOLEAN_T output_sequence, ///< should I output peptide sequence -in
  SCORER_TYPE_T output_mode  ///< the output mode -in
  )
{
  char* peptide_sequence = NULL;
  
  // print according to the output mode
  //

  SCORER_TYPE_T primary_score = output_mode;
  SCORER_TYPE_T secondary_score = SP;
  switch (output_mode) {
    case DOTP:
      // FIXME fill in once implemented
      break;
    case SP:
    case XCORR:
    case LOGP_EXP_SP:
      //case LOGP_BONF_EXP_SP:
    case LOGP_WEIBULL_SP:
    case LOGP_BONF_WEIBULL_SP:
    case DECOY_XCORR_QVALUE:
    case DECOY_PVALUE_QVALUE:
      //case LOGP_EVD_XCORR:
    case LOGP_BONF_EVD_XCORR:
    case LOGP_WEIBULL_XCORR:
    case LOGP_BONF_WEIBULL_XCORR:
    case LOGP_QVALUE_WEIBULL_XCORR:
      secondary_score = SP;
      break;
    case Q_VALUE:      
    case PERCOLATOR_SCORE:  
      primary_score = PERCOLATOR_SCORE;
      secondary_score = Q_VALUE;      
      break;
    case QRANKER_Q_VALUE:      
    case QRANKER_SCORE:  
      primary_score = QRANKER_SCORE;
      secondary_score = QRANKER_Q_VALUE;      
      break;
  }

  if(output_mode == Q_VALUE || output_mode == PERCOLATOR_SCORE){
    fprintf(file, "P\t%i\t%i\t%d\t%d\t%.9f\t%.9f\t%.9f\t", 
            get_spectrum_first_scan(match->spectrum),
            match->charge,
            match->match_rank[primary_score], 
            match->match_rank[primary_score], 
            get_peptide_peptide_mass(match->peptide), 
            match->match_scores[primary_score], 
            match->match_scores[secondary_score]);
  }
  else{
    fprintf(file, "P\t%i\t%i\t%d\t%d\t%.9f\t%.9f\t%.9f\t", 
            get_spectrum_first_scan(match->spectrum),
            match->charge,
            match->match_rank[primary_score], 
            match->match_rank[secondary_score], 
            get_peptide_peptide_mass(match->peptide), 
            match->match_scores[primary_score], 
            match->match_scores[secondary_score]);
  }
  
  // FIXME resolve spectrum_header output and above to not be coupled
 
  // should I print sequence?
  if(output_sequence){        
    peptide_sequence = get_match_sequence(match);
    fprintf(file, "%s\n", peptide_sequence);
    free(peptide_sequence);
  }
}

/**
 * \brief Print the match information in sqt format to the given file
 *
 * The main score goes in the position usually holding the xcorr.  The other
 * score goes in the position usually holding the preliminary Sp
 * score.  For searches analyzed by percolator, main and other should
 * be discriminant score and qvalue.  For p-value estimation, main and
 * other should be p-value and xcorr.
 */
void print_match_sqt(
  MATCH_T* match,             ///< the match to print -in  
  FILE* file,                 ///< output stream -out
  SCORER_TYPE_T main_score,   ///< the main score to report -in
  SCORER_TYPE_T other_score  ///< the other score to report -in
  ){

  if( match == NULL || file == NULL ){
    carp(CARP_ERROR, "Cannot print match to sqt file from null inputs");
    return;
  }
  PEPTIDE_T* peptide = get_match_peptide(match);
  // this should get the sequence from the match, not the peptide
  char* sequence = get_match_sequence_sqt(match);
  BOOLEAN_T adjust_delta_cn = FALSE;

  // NOTE (BF 12-Feb-08) This is an ugly fix to give post-percolator
  // sqt files the rank of the xcorr and sp.
  // ALSO deltaCn not serialized, so set to 0 for post-process
  SCORER_TYPE_T main_rank_type = main_score;
  SCORER_TYPE_T other_rank_type = other_score;
  if( (main_score == PERCOLATOR_SCORE || main_score == QRANKER_SCORE) && 
      (other_score == Q_VALUE || other_score == QRANKER_Q_VALUE )){
    main_rank_type = XCORR;
    other_rank_type = SP;    
    adjust_delta_cn = TRUE;
  }
  // for p-values, also give rank of xcorr and sp?
  else if( main_score == LOGP_BONF_WEIBULL_XCORR ){
    //other_rank_type = XCORR;
    main_rank_type = XCORR;
    other_score = XCORR;
  }else if( main_score == LOGP_BONF_WEIBULL_SP ){
    //other_rank_type = SP;
    other_score = SP;
  }
  // for post-analysis of p-values, both ranks from xcorr
  else if( main_score == LOGP_QVALUE_WEIBULL_XCORR ){
    main_rank_type = XCORR;
  }
  // secondary rank could always be preliminary score

  // NOTE (BF 12-Feb-08) Here is another ugly fix for post-analysis.
  // Only the fraction matched is serialized.  The number possible can
  // be calculated from the length of the sequence and the charge, but
  // the charge of the match is not serialized and I'm not sure when
  // it is set.  But I know it exists by here, so recalculate it.
  int b_y_total = get_match_b_y_ion_possible(match);
  int b_y_matched = get_match_b_y_ion_matched(match);
  
  if( b_y_total == 0 ){
    int factor = get_match_charge(match);
    if( factor == 3 ){
      factor = 2;  //there are +1 and +2 b/y ions for charge==3
    }else{
      factor = 1;
    }
    b_y_total = (get_peptide_length(peptide)-1) * 2 * factor;
    b_y_matched = (get_match_b_y_ion_fraction_matched(match)) * b_y_total;
  }

  FLOAT_T delta_cn = get_match_delta_cn(match);
  if( adjust_delta_cn == TRUE ){
    delta_cn = 0.0;
  }

  // If a p-value couldn't be calculated, print as NaN
  FLOAT_T score_main = get_match_score(match, main_score);
  if( main_score == LOGP_BONF_WEIBULL_XCORR&& score_main == P_VALUE_NA ){
    score_main = NaN();
  } 
  if (score_main == 0) { // Avoid -0
    score_main = 0.0;
  }

  // write format string with variable precision
  int precision = get_int_parameter("precision");

  // print match info
  fprintf(file, "M\t%i\t%i\t%.4f\t%.2f\t%.*g\t%.*g\t%i\t%i\t%s\tU\n",
            //fprintf(file, format,
          get_match_rank(match, main_rank_type),
          get_match_rank(match, other_rank_type),
          get_peptide_peptide_mass(peptide),
          delta_cn,
          precision,
          score_main,
          precision,
          get_match_score(match, other_score),
          b_y_matched,
          b_y_total,
          sequence
          );
  free(sequence);
  
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator = 
    new_peptide_src_iterator(peptide);
  PEPTIDE_SRC_T* peptide_src = NULL;
  char* protein_id = NULL;
  PROTEIN_T* protein = NULL;
  char* rand = "";
  if( match->null_peptide ){
    rand = "rand_";
  }
  
  while(peptide_src_iterator_has_next(peptide_src_iterator)){
    peptide_src = peptide_src_iterator_next(peptide_src_iterator);
    protein = get_peptide_src_parent_protein(peptide_src);
    protein_id = get_protein_id(protein);
    
    // print match info (locus line), add rand_ to locus name for decoys
    fprintf(file, "L\t%s%s\n", rand, protein_id);      
    free(protein_id);
  }
  
  free_peptide_src_iterator(peptide_src_iterator);
  
  return;
}

/**
 * \brief Print the match information in tab delimited format to the given file
 *
 */
void print_match_tab(
  MATCH_T* match,             ///< the match to print -in  
  FILE* file,                 ///< output stream -out
  int scan_num,             ///< starting scan number -in
  FLOAT_T spectrum_precursor_mz, ///< m/z of spectrum precursor -in
  FLOAT_T spectrum_mass,       ///< spectrum neutral mass -in
  int num_matches,            ///< num matches in spectrum -in
  int charge,                 ///< charge -in
  const BOOLEAN_T* scores_computed ///< scores_computed[TYPE] = T if match was scored for TYPE
  ){

  if( file == NULL ){ // usually b/c no decoy file to print to
    return;
  }

  if( match == NULL  ){
    carp(CARP_ERROR, 
         "Cannot print NULL match to tab delimited file.");
    return;
  }

  PEPTIDE_T* peptide = get_match_peptide(match);
  double peptide_mass = get_peptide_peptide_mass(peptide);
  // this should get the sequence from the match, not the peptide
  char* sequence = get_match_mod_sequence_str(match);
  if( sequence == NULL ){
    sequence = my_copy_string("");  // for post-search, no shuffled sequences
  }
  BOOLEAN_T adjust_delta_cn = FALSE;

  // NOTE (BF 12-Feb-08) Here is another ugly fix for post-analysis.
  // Only the fraction matched is serialized.  The number possible can
  // be calculated from the length of the sequence and the charge, but
  // the charge of the match is not serialized and I'm not sure when
  // it is set.  But I know it exists by here, so recalculate it.
  int b_y_total = get_match_b_y_ion_possible(match);
  int b_y_matched = get_match_b_y_ion_matched(match);
  
  if( b_y_total == 0 ){
    int factor = get_match_charge(match);
    if( factor == 3 ){
      factor = 2;  //there are +1 and +2 b/y ions for charge==3
    }else{
      factor = 1;
    }
    b_y_total = (get_peptide_length(peptide)-1) * 2 * factor;
    b_y_matched = (get_match_b_y_ion_fraction_matched(match)) * b_y_total;
  }

  FLOAT_T delta_cn = get_match_delta_cn(match);
  if( adjust_delta_cn == TRUE ){
    delta_cn = 0.0;
  }
  if( delta_cn == 0 ){// I hate -0, this prevents it
    delta_cn = 0.0;
  }

  BOOLEAN_T sp_scored = scores_computed[SP];
  double sp_score = get_match_score(match, SP);
  int  sp_rank = get_match_rank(match, SP);
  double xcorr_score = get_match_score(match, XCORR);
  int xcorr_rank = get_match_rank(match, XCORR);
  double log_pvalue = get_match_score(match, LOGP_BONF_WEIBULL_XCORR);
  double weibull_qvalue = get_match_score(match, LOGP_QVALUE_WEIBULL_XCORR);
  double decoy_x_qvalue = get_match_score(match, DECOY_XCORR_QVALUE);
  double decoy_p_qvalue = get_match_score(match, DECOY_PVALUE_QVALUE);
  double percolator_score = get_match_score(match, PERCOLATOR_SCORE);
  double percolator_rank = get_match_rank(match, Q_VALUE);
  double percolator_qvalue = get_match_score(match, Q_VALUE);
  double qranker_score = get_match_score(match, QRANKER_SCORE);
  double qranker_qvalue = get_match_score(match, QRANKER_Q_VALUE);
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* dig_str = digest_type_to_string(digestion);
  char *protein_ids = get_protein_ids(peptide);
  char *flanking_aas = get_flanking_aas(peptide);

  int precision = get_int_parameter("precision");
  char float_format[16];
  sprintf(float_format, "%%.%ig\t", precision);

  // Print tab delimited fields
  fprintf(file, "%d\t", scan_num);
  fprintf(file, "%d\t", charge);
  fprintf(file, "%.4f\t", spectrum_precursor_mz);
  fprintf(file, "%.4f\t", spectrum_mass);
  fprintf(file, "%.4f\t", peptide_mass);
  fprintf(file, float_format, delta_cn);
  if (sp_scored == FALSE){
    fprintf(file, "\t\t"); //score and rank
  }else{
    fprintf(file, float_format, sp_score);
    fprintf(file, "%d\t", sp_rank);
  }
  fprintf(file, float_format, xcorr_score);
  fprintf(file, "%d\t", xcorr_rank);
  if( scores_computed[LOGP_BONF_WEIBULL_XCORR] == TRUE ){ 
    // print p-value
    if (P_VALUE_NA == log_pvalue) {
      fprintf(file, "NaN\t");
    }
    else {
      fprintf(file, float_format, exp(-1 * log_pvalue));
    }
  }
  else {
    // no p-value
    fprintf(file, "\t");
  }
  if( scores_computed[LOGP_QVALUE_WEIBULL_XCORR] == TRUE ){ 
    // print q-value (Weibull est.)
    fprintf(file, float_format, exp(-1 * weibull_qvalue));
  }
  else {
    fprintf(file, "\t");
  }
  if( scores_computed[DECOY_XCORR_QVALUE]  && match->null_peptide == FALSE ){
    fprintf(file, float_format, decoy_x_qvalue);
  }
  else {
    fprintf(file, "\t");
  }
  if( scores_computed[DECOY_PVALUE_QVALUE] && match->null_peptide == FALSE){
    if (P_VALUE_NA == decoy_p_qvalue) {
      fprintf(file, "NaN\t");
    }else{
      fprintf(file, float_format, decoy_p_qvalue);
    }
  }else {
    fprintf(file, "\t");
  }
  if (scores_computed[PERCOLATOR_SCORE] == TRUE)  {
    // print percolator score
    fprintf(file, float_format, percolator_score);
    // print percolator rank
    fprintf(file, float_format, percolator_rank);
    // print q-value
    fprintf(file, float_format, percolator_qvalue);
  }
  else {
    // no percolator score, rank, or p-value
    fprintf(file, "\t\t\t");
  }
  // Output of q-ranker score and q-value will be handled here where available.
  // For now always print an empty column. 
  if (scores_computed[QRANKER_SCORE] == TRUE) {
    // print q-ranker score
    fprintf(file, float_format, qranker_score);
    // print q-value
    fprintf(file, float_format, qranker_qvalue);
  }
  else {
    fprintf(file, "\t\t");
  }
  if (sp_scored == 0 ){
    fprintf(file, "\t");
  }else{
    fprintf(file, "%d\t", b_y_matched);
  }
  fprintf(file, "%d\t", b_y_total);
  fprintf(file, "%d\t", num_matches); // Matches per spectrum
  fprintf(file, "%s\t", sequence);
  fprintf(file, "%s-%s\t", enz_str, dig_str);
  fprintf(file, "%s\t%s", protein_ids, flanking_aas);

  // if the peptide is a decoy, print the unshuffled version of the peptide
  if(match->null_peptide == TRUE){
    char* seq = get_peptide_modified_sequence(match->peptide);
    fprintf(file, "\t%s", seq);
    free(seq);
  }else{
    fprintf(file, "\t");
  }
  // End record
  fputc('\n', file);
  
  free(flanking_aas);
  free(protein_ids);
  free(sequence);
  free(enz_str);
  free(dig_str);
  
  return;
}

/**
 * shuffle the matches in the array between index start and end-1
 */
void shuffle_matches(
  MATCH_T** match_array, ///< the match array to shuffle  
  int start_index,       ///< index of first element to shuffle
  int end_index          ///< index AFTER the last element to shuffle
  ){
  if( match_array == NULL ){
    carp(CARP_ERROR, "Cannot shuffle null match array.");
    return;
  }
  //  srand(time(NULL));

  int match_idx = 0;
  for(match_idx = start_index; match_idx < end_index-1; match_idx++){
    MATCH_T* cur_match = match_array[match_idx];

    // pick a random index between match_index and end_index-1
    int rand_idx = get_random_number_interval(match_idx+1, end_index-1);

    //    fprintf(stderr, "%i values between %i and %i, rand %.4f, index %i\n",
    //            range, match_idx, end_index, rand_scaler, rand_idx);
    match_array[match_idx] = match_array[rand_idx];    
    match_array[rand_idx] = cur_match;
  }
}


/**
 * sort the match array with the corresponding compare method
 */
void qsort_match(
  MATCH_T** match_array, ///< the match array to sort -in  
  int match_total,  ///< the total number of match objects -in
  void* compare_method ///< the compare method to use -in
  )
{
  qsort(match_array, match_total, sizeof(MATCH_T*), compare_method);
}

// START need to rewrite parse_match and serialize_match to include 
// additional features

/**
 * \brief Writes the match to file in binary.
 *
 * <PEPTIDE_T: serialize peptide>
 * <float: score><int: ranking>* <--serialize for all score types
 * <SPECTRUM_T: serilize spectrum>
 * <float: b_y ion match ratio for SP>
 * <PEPTIDE_TYPE_T: the peptide type over-all peptide srcs>
 * <BOOLEAN_T: is this a null peptide?>
 *
 */
void serialize_match(
  MATCH_T* match, ///< the match to print -in
  FILE* file ///< output stream -out
  )
{
  // first serialize peptide
  serialize_peptide(match->peptide, file);
  
  // Serialize each score and rank
  int score_type_idx;
  // We don't want to change the CSM files contents so we omit q-ranker scores
  // which were added to Crux after the CSM file format had been established.
  int score_type_max = _SCORE_TYPE_NUM - 2;
  for(score_type_idx = 0; score_type_idx < score_type_max; ++score_type_idx){
    fwrite(&(match->match_scores[score_type_idx]), sizeof(FLOAT_T), 1, file);
    fwrite(&(match->match_rank[score_type_idx]), sizeof(int), 1, file);
  }
  
  // serialize spectrum in binary
  serialize_spectrum(match->spectrum, file);
  
  // b/y ion matches ratio
  fwrite(&(match->b_y_ion_fraction_matched), sizeof(FLOAT_T), 1, file);

  // serialize match peptide overall trypticity
  //fwrite(&(match->overall_type), sizeof(PEPTIDE_TYPE_T), 1, file);
  fwrite(&(match->digest), sizeof(DIGEST_T), 1, file);
  
  // serialize match is it null_peptide?
  fwrite(&(match->null_peptide), sizeof(BOOLEAN_T), 1, file);

}

/*******************************************
 * match post_process extension
 ******************************************/

/**
 * Constructs the 20 feature array that pass over to percolator registration
 * Go to top README for N,C terminus tryptic feature info.
 *\returns the feature FLOAT_T array
 */
double* get_match_percolator_features(
  MATCH_T* match, ///< the match to work -in 
  MATCH_COLLECTION_T* match_collection ///< the match collection to iterate -in
  )
{
  int feature_count = 20;
  PEPTIDE_SRC_ITERATOR_T* src_iterator = NULL;
  PEPTIDE_SRC_T* peptide_src = NULL;
  PROTEIN_T* protein = NULL;
  unsigned int protein_idx = 0;
  double* feature_array = (double*)mycalloc(feature_count, sizeof(double));
  FLOAT_T weight_diff = get_peptide_peptide_mass(match->peptide) -
    get_spectrum_neutral_mass(match->spectrum, match->charge);

  // Xcorr
  feature_array[0] = get_match_score(match, XCORR);
  // DeltCN
  feature_array[1] = match->delta_cn;
  // DeltLCN
  feature_array[2] = match->ln_delta_cn;
  // SP
  feature_array[3] = get_match_score(match, SP);
  // lnrSP
  feature_array[4] = logf(get_match_rank(match, SP));
  // dM
  feature_array[5] = weight_diff;
  // absdM
  feature_array[6] = fabsf(weight_diff);
  // Mass
  feature_array[7] = get_spectrum_neutral_mass(match->spectrum, match->charge);
  // ionFrac
  feature_array[8] = match->b_y_ion_fraction_matched;
  // lnSM
  feature_array[9] = match->ln_experiment_size;
  
  // peptide cleavage info.
  // START figure out the right way to set these features for on the fly
  // peptide generation
/*
  if(match->overall_type == TRYPTIC){
    feature_array[10] = TRUE;
    feature_array[11] = TRUE;
  }
  else if(match->overall_type == N_TRYPTIC){
    feature_array[10] = TRUE;
  }
  else if(match->overall_type == C_TRYPTIC){
    feature_array[11] = TRUE;
  }
  */
feature_array[10] = TRUE; // TODO(if perc support continues, figure out what these should be
feature_array[11] = TRUE;
  // get the missed cleave sites
  feature_array[12] = get_peptide_missed_cleavage_sites(match->peptide);
  
  // pepLen
  feature_array[13] = get_peptide_length(match->peptide);
  
  // set charge
  if(match->charge == 1){
    feature_array[14] = TRUE;
  }
  else if(match->charge == 2){
    feature_array[15] = TRUE;
  }
  else if(match->charge == 3){
    feature_array[16] = TRUE;
  }
  
  // run specific features
  if (strcmp(
        get_string_parameter_pointer("percolator-intraset-features"), "T")==0){
    carp(CARP_DETAILED_DEBUG, "Using intraset features!");
    feature_array[17] 
      = get_match_collection_hash(match_collection, match->peptide);
    
    src_iterator = new_peptide_src_iterator(match->peptide);
    // iterate overall parent proteins
    // find largest numProt and pepSite among the parent proteins
    while(peptide_src_iterator_has_next(src_iterator)){
      peptide_src = peptide_src_iterator_next(src_iterator);
      protein = get_peptide_src_parent_protein(peptide_src);
      protein_idx = get_protein_protein_idx(protein);
      
      // numProt
      if(feature_array[18] < get_match_collection_protein_counter(
                              match_collection, protein_idx)){
        feature_array[18] = get_match_collection_protein_counter(
                              match_collection, protein_idx);
      }
      
      // pepSite
      if(feature_array[19] < get_match_collection_protein_peptide_counter(
                              match_collection, protein_idx)){
        feature_array[19] = get_match_collection_protein_peptide_counter(
                              match_collection, protein_idx);      
      }
    }
  } else {
    feature_array[17] = feature_array[18] = feature_array[19] = 0.0;
  }
  
  // now check that no value is with in infinity
  int check_idx;
  for(check_idx=0; check_idx < 20; ++check_idx){
    FLOAT_T feature = feature_array[check_idx];
    if(feature <= -BILLION || feature  >= BILLION){
      carp(CARP_ERROR,
          "Percolator feature out of bounds: %d, with value %.2f. Modifying.",
           check_idx, feature);
      feature_array[check_idx] = feature <= - BILLION ? -BILLION : BILLION;
    }
  }
    
  free_peptide_src_iterator(src_iterator);
  
  return feature_array;
}

/**
 *\returns a match object that is parsed from the serialized result file
 */
MATCH_T* parse_match(
  FILE* result_file,  ///< the result file to parse PSMs -in
  DATABASE_T* database ///< the database to which the peptides are created -in
  // int num_top_match  ///< number of top PSMs serialized per spectrum -in
  )
{
  MATCH_T* match = new_match();
  carp(CARP_DETAILED_DEBUG, "New match charge is %d",get_match_charge(match));
  SPECTRUM_T* spectrum = NULL;
  PEPTIDE_T* peptide = NULL;
  
  // this is a post_process match object
  match->post_process_match = TRUE;
  int score_type_idx = 0;
  
  // parse peptide
  if((peptide = parse_peptide(result_file, database, TRUE))== NULL){
    carp(CARP_ERROR, "Failed to parse peptide");
    // FIXME should this exit or return null. I think sometimes we can get
    // no peptides, which is valid, in which case NULL makes sense.
    // maybe this should be fixed at the serialize match level however.
    return NULL;
  }
  carp(CARP_DETAILED_DEBUG, "Finished parsing match peptide.");
  // parse each score and rank of match
  // We don't want to change the CSM files contents so we omit q-ranker scores
  // which were added to Crux after the CSM file format had been established.
  int score_type_max = _SCORE_TYPE_NUM - 2;
  for(score_type_idx=0; score_type_idx < score_type_max; ++score_type_idx){
    fread(&(match->match_scores[score_type_idx]), 
      sizeof(FLOAT_T), 1, result_file);
    fread(&(match->match_rank[score_type_idx]), 
      sizeof(int), 1, result_file);
  }
  
  // parse spectrum
  if((spectrum = parse_spectrum_binary(result_file))== NULL){
    carp(CARP_ERROR, "Failed to parse binary spectrum.");
  }
  
  // spectrum specific features
  fread(&(match->b_y_ion_fraction_matched), sizeof(FLOAT_T), 1, result_file);

  // calculate the total matched from the total possible and the
  // fraction matched
  // We could do this if we had the charge. It mysteriously appears
  // later.  Calcualte this then???
  /*int total_ions = (get_peptide_length(peptide) -1) * 2;
  int matched_ions = match->b_y_ion_fraction_matched * total_ions ;
  match->b_y_ion_matched = matched_ions;
  match->b_y_ion_possible = total_ions;
  int scan = get_spectrum_first_scan(spectrum);
  carp(CARP_DETAILED_DEBUG, "For scan# %d, %d matched of %i = %.2f", 
  scan, matched_ions, total_ions, match->b_y_ion_fraction_matched);*/

  // parse match peptide overall trypticity
  //fread(&(match->overall_type), sizeof(PEPTIDE_TYPE_T), 1, result_file);
  fread(&(match->digest), sizeof(DIGEST_T), 1, result_file);
  
  // parse if match is it null_peptide?
  fread(&(match->null_peptide), sizeof(BOOLEAN_T), 1, result_file);

  // assign fields
  match->peptide_sequence = NULL;
  match->spectrum = spectrum;
  match->peptide = peptide;
  carp(CARP_DETAILED_DEBUG, "End of parse match charge is %d",
       get_match_charge(match));
  
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
char* get_match_sequence(
  MATCH_T* match ///< the match to work -in
  )
{
  // if post_process_match and has a null peptide you can't get sequence
  if(match->post_process_match && match->null_peptide){
    carp(CARP_ERROR, 
         "Cannot retrieve null peptide sequence for post_process_match");
    return NULL;
  }
  
  // if peptide sequence is cached
  // return copy of cached peptide sequence
  if(match->peptide_sequence != NULL){
    return my_copy_string(match->peptide_sequence);
  }
  
  // if not cached go generate the sequence

  // First, is this a null peptide?
  // Then must use the shuffled sequence
  if(match->null_peptide){
    // generate the shuffled peptide sequence
    if( get_boolean_parameter("reverse-sequence") == TRUE ){
      match->peptide_sequence = generate_reversed_sequence(match->peptide);
    }else{
      match->peptide_sequence = 
        generate_shuffled_sequence(match->peptide);
    }
    IF_CARP_DETAILED_DEBUG(
      char* seq = get_peptide_sequence(match->peptide);
      carp(CARP_DETAILED_DEBUG, "Shuffling transforms: %s -> %s", 
	   seq, match->peptide_sequence);
      free(seq);
    )
  }
  else{
    // just go parse it out from protein, no need to shuffle
    match->peptide_sequence = get_peptide_sequence(match->peptide);
  }
  
  return my_copy_string(match->peptide_sequence); 
  // return match->peptide_sequence;
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
char* get_match_sequence_sqt(
  MATCH_T* match ///< the match to work -in
  ){
  if( match == NULL ){
    carp(CARP_ERROR, "Cannot return sequence from NULL match");
    return NULL;
  }
  // get_match_mod_sequence (use method in case match->mod_seq == NULL) 
  MODIFIED_AA_T* mod_seq = get_match_mod_sequence(match);
  if( mod_seq == NULL ){
    return NULL;
  }
  int length = get_peptide_length(get_match_peptide(match));

  // turn it into string
  char* seq = modified_aa_string_to_string(mod_seq, length);

  // get peptide flanking residues 
  char c_term = get_peptide_c_term_flanking_aa(match->peptide);
  char n_term = get_peptide_n_term_flanking_aa(match->peptide);

  // allocate seq + 4 length array
  char* final_string = (char*)mycalloc((strlen(seq)+5), sizeof(char));

  // copy pieces in
  final_string[0] = c_term;
  final_string[1] = '.';
  strcpy(&final_string[2], seq);
  final_string[strlen(seq) + 2] = '.';
  final_string[strlen(seq) + 3] = n_term;
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
MODIFIED_AA_T* get_match_mod_sequence(
  MATCH_T* match ///< the match from which to get the sequence -in
  )
{
  if( match == NULL ){
    carp(CARP_ERROR, "Cannot get mod sequence from null match.");
  }
  // if post_process_match and has a null peptide you can't get sequence
  if(match->post_process_match && match->null_peptide){
    return NULL;
  }

  int length = get_peptide_length(get_match_peptide(match));
  
  // if peptide sequence is cached
  // return copy of cached peptide sequence
  if(match->mod_sequence != NULL){
    return copy_mod_aa_seq(match->mod_sequence, length);
  }

  // if not cached generate the sequence

  // Is this a null peptide? Then shuffle the sequence
  if(match->null_peptide){
    // generate the shuffled peptide sequence
    if( get_boolean_parameter("reverse-sequence") == TRUE){
      match->mod_sequence = generate_reversed_mod_sequence(match->peptide);
    }else{
      match->mod_sequence =
        generate_shuffled_mod_sequence(match->peptide);
    }
    IF_CARP_DETAILED_DEBUG(
      char* seq = get_peptide_sequence(match->peptide);
      char* modseq = modified_aa_string_to_string(match->mod_sequence, length);
      carp(CARP_DETAILED_DEBUG, "Shuffling transforms: %s -> %s",
	   seq, modseq );
      free(modseq);
      free(seq);
    )
  }
  else{
    // just get it from the peptide, no need to shuffle
    match->mod_sequence = get_peptide_modified_aa_sequence(match->peptide);
  }

  return copy_mod_aa_seq(match->mod_sequence, length);
}

/**
 * \brief Returns a newly allocated string of sequence including any
 * modification characters. 
 * \returns The peptide sequence of the match including modification
 * characters. 
 */
char* get_match_mod_sequence_str( MATCH_T* match ){

  // if post_process_match and has a null peptide you can't get sequence
  if(match->post_process_match && match->null_peptide){
    return NULL;
  }

  int length = get_peptide_length(get_match_peptide(match));
  
  // if sequence is cached return copy of cached peptide sequence
  if(match->mod_sequence != NULL){
    return modified_aa_string_to_string(match->mod_sequence, length);
  }

  // if not cached generate the sequence

  // Is this a null peptide? Then shuffle the sequence
  if(match->null_peptide){
    // generate the shuffled peptide sequence
    match->mod_sequence =
      generate_shuffled_mod_sequence(match->peptide);//, match->overall_type);
    IF_CARP_DETAILED_DEBUG(
      char* seq = get_peptide_sequence(match->peptide);
      char* modseq = modified_aa_string_to_string(match->mod_sequence, length);
      carp(CARP_DETAILED_DEBUG, "Shuffling transforms: %s -> %s",
	   seq, modseq );
      free(modseq);
      free(seq);
    )
  }
  else{
    // just get it from the peptide, no need to shuffle
    match->mod_sequence = get_peptide_modified_aa_sequence(match->peptide);
  }

  return modified_aa_string_to_string(match->mod_sequence, length);
}

/**
 * Must ask for score that has been computed
 *\returns the match_mode score in the match object
 */
FLOAT_T get_match_score(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  assert(match != NULL );
  return match->match_scores[match_mode];
}

/**
 * sets the match score
 */
void set_match_score(
  MATCH_T* match, ///< the match to work -out
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  FLOAT_T match_score ///< the score of the match -in
  )
{
  match->match_scores[match_mode] = match_score;
}

/**
 * Must ask for score that has been computed
 *\returns the match_mode rank in the match object
 */
int get_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode ///< the working mode (SP, XCORR) -in
  )
{
  return match->match_rank[match_mode];
}

/**
 * sets the rank of the match
 */
void set_match_rank(
  MATCH_T* match, ///< the match to work -in  
  SCORER_TYPE_T match_mode, ///< the working mode (SP, XCORR) -in
  int match_rank ///< the rank of the match -in
  )
{
  match->match_rank[match_mode] = match_rank;
}

/**
 *\returns the spectrum in the match object
 */
SPECTRUM_T* get_match_spectrum(
  MATCH_T* match ///< the match to work -in  
  )
{
  return match->spectrum;
}

/**
 * sets the match spectrum
 */
void set_match_spectrum(
  MATCH_T* match, ///< the match to work -out
  SPECTRUM_T* spectrum  ///< the working spectrum -in
  )
{

  match->spectrum = spectrum;
}

/**
 *\returns the peptide in the match object
 */
PEPTIDE_T* get_match_peptide(
  MATCH_T* match ///< the match to work -in  
  )
{
  return match->peptide;
}

/**
 * sets the match peptide, and also determines the peptide trypticity 
 * Go to top README for N,C terminus tryptic feature info.
 */
void set_match_peptide(
  MATCH_T* match, ///< the match to work -out
  PEPTIDE_T* peptide  ///< the working peptide -in
  )
{
  // first set peptide
  match->peptide = peptide;

// overall trypticity already set in peptide
//match->digest = get_peptide_digest(peptide);
  match->digest = NON_SPECIFIC_DIGEST;  // FIXME
/*
  // now set peptide overall trypticity
  PEPTIDE_SRC_ITERATOR_T* src_iterator = 
    new_peptide_src_iterator(peptide);

  PEPTIDE_SRC_T* peptide_src = NULL;

  // iterate overall parent proteins
  // determine the match overal trypticity
  // for more detail look at README at top
  while(peptide_src_iterator_has_next(src_iterator)){
    peptide_src = peptide_src_iterator_next(src_iterator);

    // now if its tryptic we are done
    if(get_peptide_src_peptide_type(peptide_src) == TRYPTIC){
      match->overall_type = TRYPTIC;
      break;
    }
    else if(get_peptide_src_peptide_type(peptide_src) == N_TRYPTIC){
      // now there're at least tryptic on both side
      // thus now we set as Tryptic and are done
      if(match->overall_type == C_TRYPTIC){
        match->overall_type = TRYPTIC;
        break;
      }
      else{
        match->overall_type = N_TRYPTIC;
      }
    }
    else if(get_peptide_src_peptide_type(peptide_src) == C_TRYPTIC){
      // now there're at least tryptic on both side
      // thus now we set as Tryptic and are done
      if(match->overall_type == N_TRYPTIC){
        match->overall_type = TRYPTIC;
        break;
      }
      else{
        match->overall_type = C_TRYPTIC;
      }
    } 
    else{
      match->overall_type = NOT_TRYPTIC;
    }
  }
  
  free_peptide_src_iterator(src_iterator);
  */
}

/**
 * sets the match if it is a null_peptide match
 */
void set_match_null_peptide(
  MATCH_T* match, ///< the match to work -out
  BOOLEAN_T is_null_peptide  ///< is the match a null peptide? -in
  )
{
  match->null_peptide = is_null_peptide;  
}

/**
 * gets the match if it is a null_peptide match
 *\returns TRUE if match is null peptide, else FALSE
 */
BOOLEAN_T get_match_null_peptide(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->null_peptide;
}

/**
 * sets the match charge
 */
void set_match_charge(
  MATCH_T* match, ///< the match to work -out
  int charge  ///< the charge of spectrum -in
  )
{
  match->charge = charge;
}

/**
 * gets the match charge
 */
int get_match_charge(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->charge;
}

/**
 * sets the match delta_cn
 */
void set_match_delta_cn(
  MATCH_T* match, ///< the match to work -out
  FLOAT_T delta_cn  ///< the delta cn value of PSM -in
  )
{
  match->delta_cn = delta_cn;
}

/**
 * gets the match delta_cn
 */
FLOAT_T get_match_delta_cn(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->delta_cn;
}

/**
 * sets the match ln_delta_cn
 */
void set_match_ln_delta_cn(
  MATCH_T* match, ///< the match to work -out
  FLOAT_T ln_delta_cn  ///< the ln delta cn value of PSM -in
  )
{
  match->ln_delta_cn = ln_delta_cn;
}

/**
 * gets the match ln_delta_cn
 */
FLOAT_T get_match_ln_delta_cn(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->ln_delta_cn;
}

/**
 * sets the match ln_experiment_size
 */
void set_match_ln_experiment_size(
  MATCH_T* match, ///< the match to work -out
  FLOAT_T ln_experiment_size ///< the ln_experiment_size value of PSM -in
  )
{
  match->ln_experiment_size = ln_experiment_size;
}

/**
 * gets the match ln_experiment_size
 */
FLOAT_T get_match_ln_experiment_size(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->ln_experiment_size;
}

/**
 * sets the match b_y_ion information
 */
void set_match_b_y_ion_info(
  MATCH_T* match, ///< the match to work -out
  SCORER_T* scorer ///< the scorer from which to extract information -in
  )
{
  match->b_y_ion_fraction_matched = 
    get_scorer_sp_b_y_ion_fraction_matched(scorer); 
  match->b_y_ion_matched = 
    get_scorer_sp_b_y_ion_matched(scorer); 
  match->b_y_ion_possible = 
    get_scorer_sp_b_y_ion_possible(scorer); 
}

/**
 * gets the match b_y_ion_fraction_matched
 */
FLOAT_T get_match_b_y_ion_fraction_matched(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->b_y_ion_fraction_matched;
}

/**
 * gets the match b_y_ion_matched
 */
int get_match_b_y_ion_matched(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->b_y_ion_matched;
}

/**
 * gets the match b_y_ion_possible
 */
int get_match_b_y_ion_possible(
  MATCH_T* match ///< the match to work -out
  )
{
  return match->b_y_ion_possible;
}

/**
 *Increments the pointer count to the match object
 */
void increment_match_pointer_count(
  MATCH_T* match ///< the match to work -in  
  )
{
  ++match->pointer_count;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

