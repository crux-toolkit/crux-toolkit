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
#include "io/carp.h"
#include "Spectrum.h"
#include "io/SpectrumCollection.h"
#include "Ion.h"
#include "IonSeries.h"
#include "util/crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "Scorer.h" 
#include "Match.h"
#include "PeptideSrc.h"
#include "ProteinIndex.h"
#include "util/modifications.h"
#include "ModifiedPeptidesIterator.h"
#include "io/MatchFileWriter.h"
#include "MatchIterator.h"
#include "io/PepXMLWriter.h"

using namespace std;

static const int _PSM_SAMPLE_SIZE = 500;
///< max number of peptides a single match collection can hold


class MatchCollection {
 friend class MatchIterator;
 protected:
  std::vector<Crux::Match*> match_; ///< array of match object
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
  vector<FLOAT_T> xcorrs_; ///< xcorrs to be used for weibull

  // The following features (post_*) are only valid when
  // post_process_collection boolean is true 
  bool post_process_collection_; 
  ///< Is this a post process match_collection?
  bool post_scored_type_set_; 
  ///< has the scored type been confirmed for the match collection,
  // set after the first match collection is extended
  Crux::Match* top_scoring_sp_; ///< the match with Sp rank == 1

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
   * Given a hash table that maps from a score to its q-value, assign
   * q-values to all of the matches in a given collection.
   */
  void assignQValues(
    const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
    SCORER_TYPE_T score_type,
    SCORER_TYPE_T derived_score_type
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
