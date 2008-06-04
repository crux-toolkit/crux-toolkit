/**
 * \file match_collection.h 
 * $Revision: 1.29.2.1 $
 * \brief A set of peptide spectrum matches for one spectrum.
 *
 * Object for given a database and a spectrum, generate all match objects
 * Creating a match collection generates all matches (searches a
 * spectrum against a database.
 */
#ifndef MATCH_COLLECTION_H
#define MATCH_COLLECTION_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>
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
#include "index.h"
#include "generate_peptides_iterator.h" 
#include "match.h"
#include "hash.h"
#include "peptide_src.h"
#include "protein_index.h"
#include "PercolatorCInterface.h"


#define _MAX_NUMBER_PEPTIDES 10000000
///< max number of peptides a single match collection can hold

// TODO (BF 1-28-08): should this be in m_c.h ?
/**
 * \returns An (empty) match_collection object.
 */
MATCH_COLLECTION_T* allocate_match_collection(void);


/**
 * \brief Creates a new match collection by searching a database
 * for matches to a spectrum. in .h
 *
 * \detail This is the main spectrum searching routine.  Allocates memory for
 * the match collection. Creates a peptide iterator for given mass
 * window. Performs preliminary scoring on all candidate
 * peptides. Performs primary scoring on the <max_rank> best-scoring
 * peptides. Estimates EVD parameters. in .h
 *
 * \returns A new match_collection object that is scored by score_type
 * and contains the top max_rank matches in .h
 */
MATCH_COLLECTION_T* new_match_collection_from_spectrum(
 SPECTRUM_T* spectrum, ///< the spectrum to match peptides in .h -in
 int charge,       ///< the charge of the spectrum -in
 int max_rank,     ///< max number of top rank matches to keep from SP -in
 SCORER_TYPE_T prelim_score, ///< the preliminary score type (SP) -in
 SCORER_TYPE_T score_type, ///< the score type (XCORR, LOGP_EXP_SP) -in
 float mass_offset,
 ///< the mass offset from neutral_mass to search for candidate peptides -in
 BOOLEAN_T null_peptide_collection,
 ///< is this match_collection a null peptide collection? -in
 INDEX_T* index,
 DATABASE_T* database
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
  MATCH_COLLECTION_T* match_collection, ///< match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  );

/**
 * \brief Sort a match_collection by the given score type, grouping
 * matches by spectrum (if multiple spectra present).
 * \returns TRUE if sort is successful, else FALSE;
 */
BOOLEAN_T spectrum_sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type to sort by -in
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
 * sets the score_type to value
 */
void set_match_collection_scored_type(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -in
  SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  BOOLEAN_T value
);

/**
 * Samples count_max matches randomly from the match_collection
 */
MATCH_COLLECTION_T* random_sample_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to sample -out
  int count_max ///< the number of matches to randomly select -in
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
 * \brief Get the number of proteins in the database associated with
 * this match collection.
 */
int get_match_collection_num_proteins(MATCH_COLLECTION_T* match_collection);

/**
 * Takes the values of match-output-folder, ms2 filename (soon to be
 * named output file), overwrite, and number-decoy-set from parameter.c 
 * and returns an array of filehandles to the newly opened files
 */
FILE** create_psm_files();

/*
 * Copied from spectrum_collection::serialize_header
 * uses values from paramter.c rather than taking as arguments
 */
void serialize_headers(FILE** file_array);

void print_matches
(MATCH_COLLECTION_T* match_collection, 
 SPECTRUM_T* spectrum, 
 BOOLEAN_T is_decoy,
 FILE* psm_file,
 FILE* sqt_file, 
 FILE* decoy_file);
/**
 * Serialize the psm features to ouput file upto 'top_match' number of 
 * top peptides among the match_collection
 *
 *
 * spectrum specific features
 * first, serialize the spectrum info of the match collection    
 * Second, iterate over matches and serialize the structs
 *
 *<int: charge state of the spectrum>
 *<int: The total match objects in the match_collection searched with the spectrum
 *<float: delta_cn>
 *<float: ln_dleta_cn>
 *<float: ln_experiment_size>
 *<BOOLEAN_T: did the score type been scored?>* <- for all score types
 *<MATCH: serialize match struct> *<--serialize match structs upto top-match # ranks
 *
 *
 *\returns TRUE, if sucessfully serializes the PSMs, else FALSE 
 */
BOOLEAN_T serialize_psm_features(
  MATCH_COLLECTION_T* match_collection, ///< working match collection -in
  FILE* output,  ///< ouput file handle -out
  int top_match, ///< number of top match to serialize -in
  SCORER_TYPE_T prelim_score, ///< the preliminary score to report -in
  SCORER_TYPE_T main_score ///<  the main score to report -in
  );

/*
 * Print the SQT file header 
 */
void print_sqt_header(FILE* outfile, char* type, int proteins);

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
 * \brief Create a match iterator to return matches from a collection
 * grouped by spectrum and sorted by given score type.
 *
 * \returns A heap-allocated match iterator.
 */
MATCH_ITERATOR_T* new_match_iterator_spectrum_sorted(
  MATCH_COLLECTION_T* match_collection, ///< match collection to iterate -in
  SCORER_TYPE_T scorer ///< the score to sort by (MATCH_SP, MATCH_XCORR) -in
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
 * Match ranks are also populated to preserve the original order of the
 * match input TRUE for preserve_order.
 *\returns TRUE, if successfully fills the scores into match object, else FALSE.
 */
BOOLEAN_T fill_result_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -out
  double* results,  ///< The result score array to fill the match objects -in
  SCORER_TYPE_T score_type,  ///< The score type of the results to fill (XCORR, Q_VALUE, ...) -in
  BOOLEAN_T preserve_order ///< preserve match order?
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
  unsigned int protein_idx ///< the protein index to return protein counter -in
  );

/**
 *\returns the match_collection protein peptide counter for the protein idx
 */
int get_match_collection_protein_peptide_counter(
  MATCH_COLLECTION_T* match_collection, ///< the working match collection -in
  unsigned int protein_idx ///< the protein index to return protein peptiide counter -in
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

/**
 * \brief Get the name of the directory the match_collection_iterator
 * is working in.
 * \returns A heap allocated string (char*) of the directory name.
 */
char* get_match_collection_iterator_directory_name(
  MATCH_COLLECTION_ITERATOR_T* iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
