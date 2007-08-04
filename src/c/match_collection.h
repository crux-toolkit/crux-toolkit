/**
 * \file match_collection.h 
 * $Revision: 1.12 $
 * \brief Object for given a database and a spectrum, generate all match objects
 */
#ifndef MATCH_COLLECTION_H
#define MATCH_COLLECTION_H

#define _MAX_NUMBER_PEPTIDES 1000000 //What to set?
#define _MIN_SP 0.010

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
 SCORER_TYPE_T score_type, ///< the score type (XCORR, LOGP_EXP_SP) -in
 float mass_offset,  ///< the mass offset from neutral_mass to search for candidate peptides -in
 BOOLEAN_T null_peptide_collection ///< is this match_collection a null peptide collection? -in
 );

/**
 * free the memory allocated match collection
 */
void free_match_collection(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  );

/**
 * sort the match collection by score_type(SP, XCORR, ... )
 *\returns TRUE, if successfully sorts the match_collection
 */
BOOLEAN_T sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to score -out
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
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
 *\returns the total match objects avaliable in current match_collection
 */
int get_match_collection_match_total(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 *\returns the total peptides searched in the experiment in match_collection
 */
int get_match_collection_experimental_size(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
  );

/**
 *\returns the top peptide count used in the logp_exp_sp in match_collection
 */
int get_match_collection_top_fit_sp(
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


/*******************************************
 * match_collection post_process extension
 ******************************************/

/**
 * create a new match collection from the serialized PSM output files
 *\returns a new match_collection object that is instantiated by the PSm output files
 */
MATCH_COLLECTION_T* new_match_collection_psm_output(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator, ///< the working match_collection_iterator -in
  SET_TYPE_T set_type  ///< what set of match collection are we creating? (TARGET, DECOY1~3) -in 
 );

/**
 * Fill the match objects score with the given the float array. 
 * The match object order must not have been altered since scoring.
 * The result array size must match the match_total count.
 *\returns TRUE, if successfully fills the scores into match object, else FALSE.
 */
BOOLEAN_T fill_result_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -out
  double* results,  ///< The result score array to fill the match objects -in
  SCORER_TYPE_T score_type  ///< The score type of the results to fill (XCORR, Q_VALUE, ...) -in
  );

/**
 * Process run specific features from all the PSMs
 */
void process_run_specific_features(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  );

/**
 *\returns the match_collection protein counter for the protein idx
 */
int get_match_collection_protein_counter(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  int protein_idx ///< the protein index to return protein counter -in
  );

/**
 *\returns the match_collection protein peptide counter for the protein idx
 */
int get_match_collection_protein_peptide_counter(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  int protein_idx ///< the protein index to return protein peptiide counter -in
  );

/**
 *\returns the match_collection hash value of PSMS for which this is the best scoring peptide
 */
int get_match_collection_hash(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  PEPTIDE_T* peptide  ///< the peptide to check hash value
  );

/******************************
 * match_collection_iterator
 ******************************/

/**
 * Create a match_collection iterator from a directory of serialized files
 * Only hadles up to one target and three decoy sets per folder
 *\returns match_collection iterator instantiated from a result folder
 */
MATCH_COLLECTION_ITERATOR_T* new_match_collection_iterator(
  char* output_file_directory, ///< the directory path where the PSM output files are located -in
  char* fasta_file ///< The name of the file (in fasta format) from which to retrieve proteins and peptides for match_collections. -in
  );

/**
 *\returns TRUE, if there's another match_collection to return, else return FALSE
 */
BOOLEAN_T match_collection_iterator_has_next(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  );

/**
 * free match_collection_iterator
 */
void free_match_collection_iterator(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  );

/**
 * returns the next match collection object and sets up fro the next iteration
 *\returns the next match collection object
 */
MATCH_COLLECTION_T* match_collection_iterator_next(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  );

/**
 *\returns the total number of match_collections to return
 */
int get_match_collection_iterator_number_collections(
  MATCH_COLLECTION_ITERATOR_T* match_collection_iterator ///< the working match_collection_iterator -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
