/**
 * \file MatchCollection.h 
 * $Revision: 1.38 $
 * \brief A set of peptide spectrum matches for one spectrum.
 *
 * Object for given a database and a spectrum, generate all match objects
 * Creating a match collection generates all matches (searches a
 * spectrum against a database.
 */
#ifndef MATCH_COLLECTION_H
#define MATCH_COLLECTION_H


#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
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
#include "Scorer.h" 
#include "Index.h"
#include "Match.h"
#include "hash.h"
#include "PeptideSrc.h"
#include "ProteinIndex.h"
#include "modifications.h"
#include "ModifiedPeptidesIterator.h"
#include "MatchFileWriter.h"
#include "MatchIterator.h"
#include "PepXMLWriter.h"

using namespace std;

static const int _PSM_SAMPLE_SIZE = 500;
static const int _MAX_NUMBER_PEPTIDES = 10000000;
///< max number of peptides a single match collection can hold


class MatchCollection {
 friend class MatchIterator;
 protected:
  Crux::Match* match_[_MAX_NUMBER_PEPTIDES]; ///< array of match object
  int match_total_;      ///< size of match array, may vary w/truncation
  int experiment_size_;  ///< total matches before any truncation
  int target_experiment_size_; ///< total target matches for same spectrum
  SpectrumZState zstate_; ///< zstate of the associated spectrum
  bool null_peptide_collection_; ///< are the peptides shuffled
  bool scored_type_[NUMBER_SCORER_TYPES]; 
                        ///< TRUE if matches have been scored by the type
  SCORER_TYPE_T last_sorted_; 
    ///< the last type by which it's been sorted ( -1 if unsorted)
  bool iterator_lock_; 
    ///< has an itterator been created? if TRUE can't manipulate matches

  bool has_distinct_matches_; ///< does the match collection have distinct matches?

  // values used for various scoring functions.
  // TODO this should be moved to match
  FLOAT_T delta_cn_; ///< the difference in top and second Xcorr scores
  FLOAT_T sp_scores_sum_; ///< for getting mean, backward compatible
  FLOAT_T sp_scores_mean_;  ///< the mean value of the scored peptides sp score
  FLOAT_T mu_;// obsolete 
  ///< EVD parameter Xcorr(characteristic value of extreme value distribution)
  FLOAT_T l_value_; // obsolete
  ///< EVD parameter Xcorr(decay constant of extreme value distribution)
  int top_fit_sp_; // obsolete
  ///< The top ranked sp scored peptides to use as EXP_SP parameter estimation
  FLOAT_T base_score_sp_; // obsolete
 ///< The lowest sp score within top_fit_sp, used as the base to rescale sp
  // Values for fitting the Weibull distribution
  FLOAT_T eta_;  ///< The eta parameter for the Weibull distribution.
  FLOAT_T beta_; ///< The beta parameter for the Weibull distribution.
  FLOAT_T shift_; ///< The location parameter for the Weibull distribution.
  FLOAT_T correlation_; ///< The correlation parameter for the Weibull distribution.
  // replace this ...
  Crux::Match* sample_matches_[_PSM_SAMPLE_SIZE];
  int num_samples_;  // the number of items in the above array
  // ...with this
  FLOAT_T xcorrs_[_MAX_NUMBER_PEPTIDES]; ///< xcorrs to be used for weibull
  int num_xcorrs_;

  // The following features (post_*) are only valid when
  // post_process_collection boolean is true 
  bool post_process_collection_; 
  ///< Is this a post process match_collection?
  bool post_scored_type_set_; 
  ///< has the scored type been confirmed for the match collection,
  // set after the first match collection is extended
  Crux::Match* top_scoring_sp_; ///< the match with Sp rank == 1

  /******* Private function declarations ***/
  int addUnscoredPeptides(
    Crux::Spectrum* spectrum, 
    SpectrumZState& charge, 
    ModifiedPeptidesIterator* peptide_iterator,
    bool is_decoy
    );

  bool scoreMatchesOneSpectrum(
    SCORER_TYPE_T score_type, 
    Crux::Spectrum* spectrum,
    int charge,
    bool store_scores
    );

  void storeNewXcorrs(
    int start_index,
    bool keep_matches
    );

  void collapseRedundantMatches();

  void consolidateMatches(
    Crux::Match** matches, 
    int start_idx, 
    int end_idx
    );

  void updateProteinCounters(
    Crux::Peptide* peptide  
    );



  /**
   * initializes a MatchCollection object
   */
  void init();


 public:
  bool exact_pval_search_;

  /**
   * \brief Creates a new match collection with no matches in it.  Sets
   * member variables from parameter.c.  The charge and null_collection
   * variables are set with the method add_matches().  Search is
   * conducted in add_matches().
   *
   * \returns A newly allocated match collection with member variables set.
   */
  MatchCollection(bool is_decoy = false);

  /**
   * create a new match collection from the serialized PSM output files
   *\returns a new match_collection object that is instantiated by the PSm output files
   */
  MatchCollection(
    MatchCollectionIterator* match_collection_iterator, ///< the working match_collection_iterator -in
    SET_TYPE_T set_type  ///< what set of match collection are we creating? (TARGET, DECOY1~3) -in 
   );

  /**
   * free the memory allocated match collection
   */
  virtual ~MatchCollection();

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
  int addMatches(
    Crux::Spectrum* spectrum,  ///< compare peptides to this spectrum
    SpectrumZState& zstate,            ///< use this charge state for spectrum
    ModifiedPeptidesIterator* peptide_iterator, ///< use these peptides
    bool is_decoy,     ///< do we shuffle the peptides
    bool store_scores, ///< true means save scores in xcorrs[]
    bool do_sp_score,  ///< true means do Sp before xcorr
    bool filter_by_sp  ///< true means keep only high sp scoring psms
  );

  /**
   * sort the match collection by score_type(SP, XCORR, ... )
   *\returns true, if successfully sorts the match_collection
   */
  void sort(
    SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
    );

  /**
   * \brief Sort a match_collection by the given score type, grouping
   * matches by spectrum (if multiple spectra are present).
   */
  void spectrumSort(
    SCORER_TYPE_T score_type ///< the score type to sort by -in
    );

  /**
   * Rank matches in a collection based on the given score type.  
   * Requires that match_collection was already scored for that score type.
   * Rank 1, means highest score
   * \returns true, if populates the match rank in the match collection
   */
  bool populateMatchRank(
    SCORER_TYPE_T score_type ///< the score type (SP, XCORR) -in
    );

  /**
   * \brief Use the matches collected from all spectra to compute FDR
   * and q_values from the ranked list of target and decoy scores.
   * Requires that matches have been scored for the given score type.
   * \returns true if q-values successfully computed, else false.
   */
  bool computeDecoyQValues();

  /**
   * match_collection get, set method
   */
  
  /**
   *\returns true, if the match collection has been scored by score_type
   */
  bool getScoredType(
    SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  );

  /**
   * sets the score_type to value
   */
  void setScoredType(
    SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
    bool value
  );

  void getCustomScoreNames(
    std::vector<std::string>& custom_score_names
  );

  void preparePostProcess();

  bool extendTabDelimited(
    Database* database, ///< the database holding the peptides -in
    MatchFileReader& result_file,   ///< the result file to parse PSMs -in
    Database* decoy_database = NULL ///< optional database with decoy peptides
    );

  /**
   * Samples count_max matches randomly from the match_collection
   */
  MatchCollection* randomSample(
    int count_max ///< the number of matches to randomly select -in
    );

  /**
   * This function is a transformation of the partial derivatives of
   * the log likelihood of the data given an extreme value distribution
   * with location parameter mu and scale parameter 1/L. The transformation 
   * has eliminated the explicit dependence on the location parameter, mu, 
   * leaving only the scale parameter, 1/L.
   *
   * The zero crossing of this function will correspond to the maximum of the 
   * log likelihood for the data.
   *
   * See equations 10 and 11 of "Maximum Likelihood fitting of extreme value 
   * distributions".
   *
   * The parameter values contains a list of the data values.
   * The parameter L is the reciprocal of the scale parameters.
   *
   *\returns the final exponential values of the score and sets the value of the function and its derivative.
   */
  void constraintFunction(
    SCORER_TYPE_T score_type, ///< score_type to estimate EVD distribution -in
    FLOAT_T l_value,  ///< L value -in
    FLOAT_T* function,  ///< the output function value -out
    FLOAT_T* derivative,  ///< the output derivative value -out
    FLOAT_T* exponential_sum ///< the final exponential array sum -out
    );

  /**
   *\returns true, if there is a  match_iterators instantiated by match collection 
   */
  bool getIteratorLock();

  /**
   *\returns the total match objects avaliable in current match_collection
   */
  int getMatchTotal();

  /**
   * Sets the total peptides searched in the experiment in match_collection
   */
  void setExperimentSize(int size);

  /**
   *\returns the total peptides searched in the experiment in match_collection
   */
  int getExperimentSize();

  /**
   * Sets the total number of target peptides searched for this
   * spectrum.  Only to be used by decoy match collections.
   */
  void setTargetExperimentSize(int numMatches);

  /**
   * \returns the number of target matches that this spectrum had.
   * Different than getExperimentSize() for decoy match collections.
   */
  int getTargetExperimentSize();

  /**
   * Set the filepath for all matches in the collection
   * \returns the associated file idx
   */
  int setFilePath(
    const std::string& file_path  ///< File path to set
  );

  /**
   * \returns true if the match_collection only contains decoy matches,
   * else (all target or mixed) returns false.
   */
  bool isDecoy();

  /**
   * Keeps track of the top scoring PSM by Sp.  To be run after
   * re-ranking by Sp.
   */
  void saveTopSpMatch();

  /**
   *\returns the top peptide count used in the logp_exp_sp in match_collection
   */
  int getTopFitSp();

  /**
   *\returns the charge of the spectrum that the match collection was created
   */
  int getCharge();

  bool calculateDeltaCn();

  /**
   * Must have been scored by Xcorr, returns error if not scored by Xcorr
   *\returns the delta cn value(difference in top and second ranked Xcorr values)
   */
  FLOAT_T getDeltaCn();

  /**
   * \brief Transfer the Weibull distribution parameters, including the
   * correlation from one match_collection to another.  No check to see
   * that the parameters have been estimated.
   */
  static void transferWeibull(
    MatchCollection* from_collection,
    MatchCollection* to_collection
    );

  /**
   * \brief Remove matches from the collection so that only those of
   * rank 'max_rank' and higher for the given score type remain.
   */
  void truncate(
    int max_rank,     
    SCORER_TYPE_T score_type 
    );

  /**
   * \brief Put all the matches from the source match collection in the
   * destination. Only copies the pointers of the matches so use with
   * caution. 
   * \returns The number of matches added.
   */
  static int merge(
    MatchCollection* source, ///< matches to be moved
    MatchCollection* destination ///< move matches to here
  );

  /**
   * \brief Add a single match to a collection.
   * Only puts a copy of the pointer to the match in the
   * match_collection, does not allocate a new match.
   */
  bool addMatch(
    Crux::Match* match                        ///< add this match
  );

  /**
   * \brief Print the given match collection for one spectrum to all
   * appropriate files. 
   */
  void print(
    Crux::Spectrum* spectrum, 
    bool is_decoy,
    FILE* psm_file,
    FILE* sqt_file, 
    FILE* decoy_file,
    FILE* tab_file, 
    FILE* decoy_tab_file
  );

  /**
   * \brief Print the given match collection for several spectra to all
   * appropriate files.  Takes the spectrum information from the matches
   * in the collection.
   */
  void printMultiSpectra(
    MatchFileWriter* tab_file, 
    MatchFileWriter* decoy_tab_file
    );

  /**
   * \brief Print the given match collection for several spectra to
   * xml files only. Takes the spectrum information from the
   * matches in the collection. At least for now, prints all matches in
   * the collection rather than limiting by top-match parameter. 
   */
  void printMultiSpectraXml(
    PepXMLWriter* output
    );

  /*
   * Print the XML file header
   */ 
  static void printXmlHeader(FILE* outfile);

  /*
   * Print the SQT file header 
   */
  static void printSqtHeader(
    FILE* outfile, 
    const char* type, 
    string database,
    int num_proteins,
    bool exact_pval_search_ = false
    );

  /*
   * Print the tab delimited file header 
   */
  static void printTabHeader(
    FILE* outfile
    );

  /*
   * Print the XML file footer
   */
  static void printXmlFooter(
    FILE* outfile
    );

  /**
   * Print the psm features to output file up to 'top_match' number of
   * top peptides among the match_collection in xml file format
   * returns true, if sucessfully print xml format of the PSMs, else false
   */
  bool printXml(
    PepXMLWriter* output,
    int top_match,
    Crux::Spectrum* spectrum,
    SCORER_TYPE_T main_score
    );

  /**
   * Print the psm features to output file upto 'top_match' number of 
   * top peptides among the match_collection in sqt file format
   *\returns true, if sucessfully print sqt format of the PSMs, else false 
   */
  bool printSqt(
    FILE* output, ///< the output file -out
    int top_match, ///< the top matches to output -in
    Crux::Spectrum* spectrum ///< the spectrum to print sqt -in
    );

  /**
   * Print the psm features to output file upto 'top_match' number of 
   * top peptides among the match_collection in tab delimited file format
   *\returns true, if sucessfully print sqt format of the PSMs, else false 
   */
  bool printTabDelimited(
    MatchFileWriter* output, ///< the output file -out
    int top_match, ///< the top matches to output -in
    Crux::Spectrum* spectrum, ///< the spectrum to print sqt -in
    SCORER_TYPE_T main_score  ///< the main score to report -in
    );

  /**
   * Retrieve the calibration parameter eta.
   */
  FLOAT_T getCalibrationEta();

  /**
   * Retrieve the calibration parameter beta.
   */
  FLOAT_T getCalibrationBeta();

  /**
   * Retrieve the calibration parameter shift.
   */
  FLOAT_T getCalibrationShift();

  /**
   * Retrieve the calibration parameter correlation.
   */
  FLOAT_T getCalibrationCorr();

  bool getHasDistinctMatches();
  void setHasDistinctMatches(bool distinct_matches);
  /**
   * Print the calibration parameters eta, beta, shift and correlation
   * with tabs between.
   */
  void printCalibrationParameters(
    FILE* output ///< The output file -in
    );

  /**
   * \brief Check that a match collection has a sufficient number of
   * matches for estimating Weibull parameters.
   * \returns true if their are enough xcorrs for estimating Weibull
   * parameters or false if not.
   */
  bool hasEnoughWeibullPoints();

  bool estimateWeibullParameters(
    SCORER_TYPE_T score_type,
    int sample_count, 
    Crux::Spectrum* spectrum,
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
  bool estimateWeibullParametersFromXcorrs(
    Crux::Spectrum* spectrum,
    int charge
    );

  bool computePValues(
    FILE* output_pvalue_file
    );

  /**
   * Try setting the match collection's zstate.  Successful if the
   * current charge is 0 (i.e. hasn't yet been set) or if the current
   * charge is the same as the given value.  Otherwise, returns false
   *
   * \returns true if the match_collection's zstate was changed.
   */
  bool setZState(
    SpectrumZState& zstate ///< new zstate
  );

  /**
   * Search the given database or index using shuffled peptides and the
   * spectrum/charge in the target psm match collection.  Add those
   * scores to the target psm match collection for use in weibull
   * parameter estimation but do not save the matches
   */
  void addDecoyScores(
    Crux::Spectrum* spectrum, ///< search this spectrum
    SpectrumZState& zstate, ///< search spectrum at this charge state
    ModifiedPeptidesIterator* peptides ///< use these peptides to search
  );

  /**
   * Extract a given type of score into an array.  The array is
   * allocated here and must be freed by the caller.
   */
  FLOAT_T* extractScores(
    SCORER_TYPE_T       score_type ///< Type of score to extract.
  );

  /**
   * Given a hash table that maps from a score to its q-value, assign
   * q-values to all of the matches in a given collection.
   */
  void assignQValues(
    const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
    SCORER_TYPE_T score_type
    );

  /**
   * Given a hash table that maps from a score to its PEP, assign
   * PEPs to all of the matches in a given collection.
   */
  void assignPEPs(
    const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
    SCORER_TYPE_T score_type
    );

  /*******************************************
   * match_collection post_process extension
   ******************************************/
  /**
   * Fill the match objects score with the given the FLOAT_T array. 
   * The match object order must not have been altered since scoring.
   * The result array size must match the match_total count.
   * Match ranks are also populated to preserve the original order of the
   * match input true for preserve_order.
   *\returns true, if successfully fills the scores into match object, else false.
   */
  void fillResult(
    double* results,  ///< The result score array to fill the match objects -in
    SCORER_TYPE_T score_type,  ///< The score type of the results to fill -in
    bool preserve_order ///< preserve match order?
    );

  /**
   * Process run specific features from all the PSMs
   */
  void processRunSpecificFeatures();

// cheater functions for testing
  void forceScoredBy(SCORER_TYPE_T type);

  bool addMatchToPostMatchCollection(
    Crux::Match* match 
    );



};

/**
 * Read files in the directory and return the names of target or
 * decoy files to use for post-search commands.
 * \returns Vector parameter filled with names of target or decoy
 * files.
 */
void get_target_decoy_filenames(vector<string>& target_decoy_names,
                                DIR* directory,
                                SET_TYPE_T type);


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
