/**
 * \file match_collection.h 
 * $Revision: 1.38 $
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
#include <map>
#include <time.h>
#include "carp.h"
#include "parse_arguments.h"
#include "Spectrum.h"
#include "SpectrumCollection.h"
#include "Ion.h"
#include "IonSeries.h"
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
#include "modifications.h"
#include "modified_peptides_iterator.h"
#include "MatchFileWriter.h"

using namespace std;

static const int _PSM_SAMPLE_SIZE = 500;
static const int _MAX_NUMBER_PEPTIDES = 10000000;
///< max number of peptides a single match collection can hold

// TODO (BF 1-28-08): should this be in m_c.h ?
/**
 * \returns An (empty) match_collection object.
 */
MATCH_COLLECTION_T* allocate_match_collection(void);

/**
 * \brief Creates a new match collection with no matches in it.  Sets
 * member variables from parameter.c.  The charge and null_collection
 * variables are set with the method add_matches().  Search is
 * conducted in add_matches().
 *
 * \returns A newly allocated match collection with member variables set.
 */
MATCH_COLLECTION_T* new_empty_match_collection(BOOLEAN_T is_decoy);

/**
 * free the memory allocated match collection
 */
void free_match_collection(
  MATCH_COLLECTION_T* match_collection ///< the match collection to free -out
  );

/**
 * \brief The main search function.  All peptides in the peptide
 * iterator are compared to the spectrum and the resulting score(s)
 * are stored in a match.  All matches are stored in the
 * match_collection.  Can be called on an empty match_collection or
 * one already containing matches.  No checks to confirm that the same
 * spectrum is being searched in subsiquent calls.
 *
 * First, the prelimiary score (as in parameter.c) is used to compare
 * peptides and spectrum.  These results are then sorted and the final
 * score (as in parameter.c) is calculated on the top-match
 * (parameter.c) top matches as ranked by the preliminary score.  No
 * matches are deleted after ranking.
 *
 * When called on a match collection already containing matches, the
 * preliminary score is calculated for all new peptides.  All matches
 * (from this peptide iterator and previous) are sorted by prelim
 * score and only the top-match matches are scored for the final
 * score.  Previously scored matches are not scored twice.
 *
 * \returns The number of matches added.
 */
int add_matches(
  MATCH_COLLECTION_T* match_collection,///< add matches to this
  Spectrum* spectrum,  ///< compare peptides to this spectrum
  SpectrumZState& zstate,            ///< use this charge state for spectrum
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator, ///< use these peptides
  BOOLEAN_T is_decoy,     ///< do we shuffle the peptides
  BOOLEAN_T store_scores, ///< TRUE means save scores in xcorrs[]
  BOOLEAN_T do_sp_score,  ///< TRUE means do Sp before xcorr
  BOOLEAN_T filter_by_sp  ///< TRUE means keep only high sp scoring psms
);

/**
 * sort the match collection by score_type(SP, XCORR, ... )
 *\returns TRUE, if successfully sorts the match_collection
 */
void sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
  );

/**
 * \brief Sort a match_collection by the given score type, grouping
 * matches by spectrum (if multiple spectra are present).
 */
void spectrum_sort_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< match collection to sort -out
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  );

/**
 * Rank matches in a collection based on the given score type.  
 * Requires that match_collection was already scored for that score type.
 * Rank 1, means highest score
 * \returns TRUE, if populates the match rank in the match collection
 */
BOOLEAN_T populate_match_rank_match_collection(
 MATCH_COLLECTION_T* match_collection, ///< the match collection to populate match rank -out
 SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
 );

/**
 * \brief Use the matches collected from all spectra to compute FDR
 * and q_values from the ranked list of target and decoy scores.
 * Requires that matches have been scored for the given score type.
 * \returns TRUE if q-values successfully computed, else FALSE.
 */
BOOLEAN_T compute_decoy_q_values(MATCH_COLLECTION_T* match_collection);

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
 * \returns TRUE if the match_collection only contains decoy matches,
 * else (all target or mixed) returns FALSE.
 */
BOOLEAN_T get_match_collection_is_decoy(
  MATCH_COLLECTION_T* match_collection
);

/**
 * Keeps track of the top scoring PSM by Sp.  To be run after
 * re-ranking by Sp.
 */
void save_top_sp_match(MATCH_COLLECTION_T* matches);

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
FLOAT_T get_match_collection_delta_cn(
  MATCH_COLLECTION_T* match_collection ///< working match collection -in
);

/**
 * \brief Get the number of proteins in the database associated with
 * this match collection.
 */
int get_match_collection_num_proteins(MATCH_COLLECTION_T* match_collection);

/**
 * \brief Transfer the Weibull distribution parameters, including the
 * correlation from one match_collection to another.  No check to see
 * that the parameters have been estimated.
 */
void transfer_match_collection_weibull(
  MATCH_COLLECTION_T* from_collection,
  MATCH_COLLECTION_T* to_collection
  );

/**
 * \brief Remove matches from the collection so that only those of
 * rank 'max_rank' and higher for the given score type remain.
 */
void truncate_match_collection(
  MATCH_COLLECTION_T* match_collection, 
  int max_rank,     
  SCORER_TYPE_T score_type 
  );

/**
 * \brief Put all the matches from the source match collection in the
 * destination. Only copies the pointers of the matches so use with
 * caution. 
 * \returns The number of matches added.
 */
int merge_match_collections(
  MATCH_COLLECTION_T* source, ///< matches to be moved
  MATCH_COLLECTION_T* destination ///< move matches to here
);

/**
 * \brief Add a single match to a collection.
 * Only puts a copy of the pointer to the match in the
 * match_collection, does not allocate a new match.
 */
BOOLEAN_T add_match_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< add to here
  MATCH_T* match                        ///< add this match
);

/**
 * \brief Print the given match collection for one spectrum to all
 * appropriate files. 
 */
void print_matches
(MATCH_COLLECTION_T* match_collection, 
 Spectrum* spectrum, 
 BOOLEAN_T is_decoy,
 FILE* psm_file,
 FILE* sqt_file, 
 FILE* decoy_file,
 FILE* tab_file, 
 FILE* decoy_tab_file);

/**
 * \brief Print the given match collection for several spectra to all
 * appropriate files.  Takes the spectrum information from the matches
 * in the collection.
 */
void print_matches_multi_spectra
(MATCH_COLLECTION_T* match_collection, 
 MatchFileWriter* tab_file, 
 MatchFileWriter* decoy_tab_file);


/**
 * \brief Print the given match collection for several spectra to
 * xml files only. Takes the spectrum information from the
 * matches in the collection. At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter. 
 */
void print_matches_multi_spectra_xml(
 MATCH_COLLECTION_T* match_collection,
 FILE* output);


/*
 * Print the XML file header
 */ 
void print_xml_header(FILE* outfile);



/*
 * Print the SQT file header 
 */
void print_sqt_header(FILE* outfile, 
                      const char* type, 
                      int proteins, 
                      BOOLEAN_T is_for_match_analysis);

/*
 * Print the tab delimited file header 
 */
void print_tab_header(FILE* outfile);


/*
 * Print the XML file footer
 */
void print_xml_footer(FILE* outfile);

/**
 * Print the psm features to output file upto 'top_match' number of
 * top peptides among the match_collection in xml file format
 * returns TRUE, if sucessfully print xml format of the PSMs, else FALSE
 */


BOOLEAN_T print_match_collection_xml(
  FILE* output,
  int top_match,
  MATCH_COLLECTION_T* match_collection,
  Spectrum* spectrum,
  SCORER_TYPE_T main_score,
  int index
  );


/**
 * Print the psm features to output file upto 'top_match' number of 
 * top peptides among the match_collection in sqt file format
 *\returns TRUE, if sucessfully print sqt format of the PSMs, else FALSE 
 */
BOOLEAN_T print_match_collection_sqt(
  FILE* output, ///< the output file -out
  int top_match, ///< the top matches to output -in
  MATCH_COLLECTION_T* match_collection, ///< the match_collection to print -in
  Spectrum* spectrum ///< the spectrum to print sqt -in
  );

/**
 * Print the psm features to output file upto 'top_match' number of 
 * top peptides among the match_collection in tab delimited file format
 *\returns TRUE, if sucessfully print sqt format of the PSMs, else FALSE 
 */
BOOLEAN_T print_match_collection_tab_delimited(
  MatchFileWriter* output, ///< the output file -out
  int top_match, ///< the top matches to output -in
  MATCH_COLLECTION_T* match_collection, ///< the match_collection to print sqt -in
  Spectrum* spectrum, ///< the spectrum to print sqt -in
  SCORER_TYPE_T main_score  ///< the main score to report -in
  );

/**
 * Retrieve the calibration parameter eta.
 */
FLOAT_T get_calibration_eta (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
);

/**
 * Retrieve the calibration parameter beta.
 */
FLOAT_T get_calibration_beta (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
);

/**
 * Retrieve the calibration parameter shift.
 */
FLOAT_T get_calibration_shift (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
);

/**
 * Retrieve the calibration parameter correlation.
 */
FLOAT_T get_calibration_corr (
  MATCH_COLLECTION_T* my_collection ///< The collection -in
);

/**
 * Print the calibration parameters eta, beta, shift and correlation
 * with tabs between.
 */
void print_calibration_parameters(
  MATCH_COLLECTION_T* my_collection, ///< The collection -in
  FILE* output ///< The output file -in
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
 * Fill the match objects score with the given the FLOAT_T array. 
 * The match object order must not have been altered since scoring.
 * The result array size must match the match_total count.
 * Match ranks are also populated to preserve the original order of the
 * match input TRUE for preserve_order.
 *\returns TRUE, if successfully fills the scores into match object, else FALSE.
 */
void fill_result_to_match_collection(
  MATCH_COLLECTION_T* match_collection, ///< the match collection to iterate -out
  double* results,  ///< The result score array to fill the match objects -in
  SCORER_TYPE_T score_type,  ///< The score type of the results to fill -in
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
  const char* output_file_directory, ///< the directory path where the PSM output files are located -in
  const char* fasta_file, ///< The name of the file (in fasta format) from which to retrieve proteins and peptides for match_collections. -in
  int* decoy_count
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

/**
 * \brief Check that a match collection has a sufficient number of
 * matches for estimating Weibull parameters.
 * \returns TRUE if their are enough xcorrs for estimating Weibull
 * parameters or FALSE if not.
 */
BOOLEAN_T has_enough_weibull_points(
  MATCH_COLLECTION_T* match_collection
);

BOOLEAN_T estimate_weibull_parameters(
  MATCH_COLLECTION_T* match_collection, 
  SCORER_TYPE_T score_type,
  int sample_count, 
  Spectrum* spectrum,
  int charge
  );

/**
 * \brief Use the xcorrs saved in the match_collection to estimate the
 * weibull parameters to be used for computing p-values. 
 *
 * Requires that main score be XCORR, but with relativly few changes
 * other scores could be accomodated.
 * Implementation of Weibull distribution parameter estimation from 
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 */
BOOLEAN_T estimate_weibull_parameters_from_xcorrs(
  MATCH_COLLECTION_T* match_collection, 
  Spectrum* spectrum,
  int charge
  );

BOOLEAN_T compute_p_values(
  MATCH_COLLECTION_T* match_collection,
  FILE* output_pvalue_file
);

/**
 * Try setting the match collection's charge.  Successful if the
 * current charge is 0 (i.e. hasn't yet been set) or if the current
 * charge is the same as the given value.  Otherwise, returns false
 *
 * \returns TRUE if the match_collection's charge state was changed.
 */
/*
BOOLEAN_T set_match_collection_charge(
  MATCH_COLLECTION_T* match_collection,  ///< match collection to change
  int charge ///< new charge state
);
*/
BOOLEAN_T set_match_collection_zstate(
  MATCH_COLLECTION_T* match_collection, ///< match collection to change
  SpectrumZState& zstate ///< new zstate
);

/**
 * Search the given database or index using shuffled peptides and the
 * spectrum/charge in the target psm match collection.  Add those
 * scores to the target psm match collection for use in weibull
 * parameter estimation but do not save the matches
 */
void add_decoy_scores_match_collection(
  MATCH_COLLECTION_T* target_matches, ///< add scores to this collection
  Spectrum* spectrum, ///< search this spectrum
  int charge, ///< search spectrum at this charge state
  MODIFIED_PEPTIDES_ITERATOR_T* peptides ///< use these peptides to search
);


/**
 * Extract a given type of score into an array.  The array is
 * allocated here and must be freed by the caller.
 */
FLOAT_T* extract_scores_match_collection(
  SCORER_TYPE_T       score_type, ///< Type of score to extract.
  MATCH_COLLECTION_T* all_matches ///< add scores to this collection
);

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
void assign_match_collection_qvalues(
  const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
  SCORER_TYPE_T score_type,
  MATCH_COLLECTION_T* all_matches
);

const vector<bool>& get_match_collection_iterator_cols_in_file(
  MATCH_COLLECTION_ITERATOR_T* match_collection);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
