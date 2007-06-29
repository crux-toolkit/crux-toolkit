/**
 * \file match_collection.h 
 * $Revision: 1.6 $
 * \brief Object for given a database and a spectrum, generate all match objects
 */
#ifndef MATCH_COLLECTION_H
#define MATCH_COLLECTION_H

#define _MAX_NUMBER_PEPTIDES 1000000 //What to set?

/**
 * \returns An (empty) match_collection object.
 */
MATCH_COLLECTION_T* allocate_match_collection(void);

/**
 * create a new match collection from spectrum
 * creates a peptide iterator for given mass window
 * return the top max_rank matches, first scored by prelim_score(SP), then by score_type(XCORR, LOGP_EXP_SP)
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in
 SCORER_TYPE_T prelim_score, ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type ///< the score type (XCORR, LOGP_EXP_SP) -in
 );

/**
 * create a new match collection from spectrum
 * return the top max_rank matches, first scored by prelim_score(SP), then by score_type(XCORR, LOGP_EXP_SP);
 * uses a provided peptide iterator, MUST be a mutable iterator
 * Sets the iterator before useage.
 *\returns a new match_collection object that is scored by score_type and contains the top max_rank matches
 */
MATCH_COLLECTION_T* new_match_collection_spectrum_with_peptide_iterator(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in 
 SCORER_TYPE_T prelim_score, ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type ///< the score type (XCORR, LOGP_EXP_SP) -in
 //GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator ///< peptide iteartor to use, must set it first before use
 );

/**
 * free the memory allocated match collection
 */
void free_match_collection(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  );

/**
 * match_collection get, set method
 */

/**
 *\returns TRUE, if the match collection has been scored by score_type
 */
BOOLEAN_T get_match_collection_scored_type(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -in
  SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
);


/**
 *\returns TRUE, if there is a  match_iterators instantiated by match collection 
 */
BOOLEAN_T get_match_collection_iterator_lock(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 *\returns the total match objects in match_collection
 */
int get_match_collection_match_total(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 *\returns the charge of the spectrum that the match collection was created
 */
int get_match_collection_charge(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 * Must have been scored by Xcorr, returns error if not scored by Xcorr
 *\returns the delta cn value(difference in top and second ranked Xcorr values)
 */
float get_match_collection_delta_cn(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 * Serialize the psm features to ouput file upto 'top_match' number of 
 * top peptides among the match_collection
 *\returns TRUE, if sucessfully serializes the PSMs, else FALSE 
 */
BOOLEAN_T serialize_psm_features(
  MATCH_COLLECTION_T* match_collection, ///< working match collection -in
  SPECTRUM_T* spectrum, ///< the working spectrum -in
  FILE* output,  ///< ouput file handle -out
  int top_match, ///< number of top match to serialize -in
  SCORER_TYPE_T prelim_score, ///< the preliminary score to report -in
  SCORER_TYPE_T main_score ///<  the main score to report -in
  );

/**
 * Print the psm features to output file upto 'top_match' number of 
 * top peptides among the match_collection in sqt file format
 *\returns TRUE, if sucessfully print sqt format of the PSMs, else FALSE 
 */
BOOLEAN_T print_match_collection_sqt(
  FILE* output, ///< the output file -out
  int top_match, ///< the top matches to output -in
  int charge, ///< the charge of the of spectrum -in
  MATCH_COLLECTION_T* match_collection, ///< the match_collection to print sqt -in
  SPECTRUM_T* spectrum, ///< the spectrum to print sqt -in
  SCORER_TYPE_T prelim_score, ///< the preliminary score to report -in
  SCORER_TYPE_T main_score  ///< the main score to report -in
  );

/**
 * match_iterator routines!
 */

/**
 * create a new memory allocated match iterator
 * creates a new the generate_peptides_iterator inside the match_iterator
 *\returns a new memory allocated match iterator
 */
MATCH_ITERATOR_T* new_match_iterator(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -in
  SCORER_TYPE_T match_mode, ///< the mode to set (MATCH_SP, MATCH_XCORR) -in
  BOOLEAN_T sort_match  ///< should I return the match in sorted order?
  );

/**
 * Does the match_iterator have another match struct to return?
 * MUST set the iterator to correct mode before calling this method
 *\returns TRUE, if match iter has a next match, else False
 */
BOOLEAN_T match_iterator_has_next(
  MATCH_ITERATOR_T* match_iterator ///< the working  match iterator -in
  );

/**
 * return the next match struct!
 * MUST set the iterator to correct mode before initialially calling this method
 *\returns the match in decreasing score order for the match_mode(SCORER_TYPE_T)
 */
MATCH_T* match_iterator_next(
  MATCH_ITERATOR_T* match_iterator ///< the working match iterator -in
  );

/**
 * free the memory allocated iterator
 */
void free_match_iterator(
  MATCH_ITERATOR_T* match_iterator ///< the match iterator to free
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
