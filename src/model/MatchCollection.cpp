/*********************************************************************//**
 * \file MatchCollection.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 11/27 2006
 * \brief A set of peptide spectrum matches for one spectrum.
 *
 * Methods for creating and manipulating match_collections.   
 * Creating a match collection generates all matches (searches a
 * spectrum against a database.
 ****************************************************************************/
#include "MatchCollection.h"
#include "MatchCollectionIterator.h"
#include "MatchIterator.h"
#include <string>
#include "io/MatchFileReader.h"
#include "io/SQTReader.h"
#include "util/WinCrux.h"

using namespace std;
using namespace Crux;

/**
 * \returns An (empty) match_collection object.
 */
void MatchCollection::init() {

  match_.reserve(10 * MILLION);
  experiment_size_ = 0;
  target_experiment_size_ = 0;
  zstate_ = SpectrumZState();
  null_peptide_collection_ = false;

  // loop over to set all score type to false
  for(int score_type_idx = 0 ; 
    score_type_idx < NUMBER_SCORER_TYPES ; 
    ++score_type_idx){

    scored_type_[score_type_idx] = false;
  }

  // set last score to -1, thus nothing has been done yet
  last_sorted_ = (SCORER_TYPE_T)-1;
  iterator_lock_ = false;
  delta_cn_ = 0;
  sp_scores_sum_ = 0;
  sp_scores_mean_ = 0;
  mu_ = 0;
  l_value_ = 0;
  top_fit_sp_ = 0;
  base_score_sp_ = 0;
  eta_ = 0;
  beta_ = 0;
  shift_ = 0;
  correlation_ = 0;
  num_samples_ = 0;
  xcorrs_.reserve(10 * MILLION);

  post_process_collection_ = false;
  post_scored_type_set_ = false;
  top_scoring_sp_ = NULL;
  exact_pval_search_ = false;
  has_distinct_matches_ = false;
}

/**
 * /brief Free the memory allocated for a match collection
 * Deep free; each match is freed which, in turn, frees each spectrum
 * and peptide. 
 */
MatchCollection::~MatchCollection() {

  // decrement the pointer count in each match object
  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    Match::freeMatch(*i);
  }

  // and free the sample matches
  while(num_samples_ > 0){
    --num_samples_;
    Match::freeMatch(sample_matches_[num_samples_]);
    sample_matches_[num_samples_] = NULL;
  }

  if(top_scoring_sp_){
    Match::freeMatch(top_scoring_sp_);
  }
}

/**
 * \brief Creates a new match collection with no matches in it.  Sets
 * the member variable indicating if the matches are to real peptides
 * or to decoy (shuffled) peptides. Other member variables are set to
 * default values.  The method add_matches() can be used to search a
 * spectrum and store the matches in this collection.  
 *
 * \returns A newly allocated match collection
 */
MatchCollection::MatchCollection(
  bool is_decoy
  ){

  init();
  null_peptide_collection_ = is_decoy;
}

void MatchCollection::preparePostProcess() {

  // prepare the match_collection
  init();
  // set this as a post_process match collection
  post_process_collection_ = true;
}


/**
 * \brief Creates a new match_collection from the match collection
 * iterator. 
 *
 * Used in the post_processing extension.  Also used by
 * setup_match_collection_iterator which is called by next to find,
 * open, and parse the next psm file(s) to process.  If there are
 * multiple target psm files, it reads in all of them when set_type is
 * 0 and puts them all into one match_collection.  If the fileroot
 * parameter is non-null, only reads files with that prefix.
 *\returns A heap allocated match_collection.
 */
MatchCollection::MatchCollection(
  MatchCollectionIterator* match_collection_iterator, 
    ///< the working match_collection_iterator -in
  SET_TYPE_T set_type  
    ///< what set of match collection are we creating? (TARGET, DECOY1~3) -in 
  )
{ 

  Database* database = match_collection_iterator->getDatabase();
  Database* decoy_database = match_collection_iterator->getDecoyDatabase();

  preparePostProcess();

  // get the list of files to open
  vector<string> file_names;
  get_target_decoy_filenames(file_names, 
                             match_collection_iterator->getWorkingDirectory(),
                             set_type);

  // open each file and add psms to match collection
  for(int file_idx = 0; file_idx < (int)file_names.size(); file_idx++){
    char* full_filename = 
      get_full_filename(match_collection_iterator->getDirectoryName(),
                        file_names[file_idx].c_str());
    MatchFileReader delimited_result_file(full_filename);
    carp(CARP_DEBUG, "Creating new match collection from '%s' file.",
         full_filename);
    free(full_filename);

    extendTabDelimited(database, delimited_result_file, decoy_database);

    // for the first target file, set headers based on input files
    if( set_type == SET_TARGET && file_idx == 0 ){
      delimited_result_file.getMatchColumnsPresent(
                                  match_collection_iterator->getColsInFile()); 
    }
  } // next file
}

/**
 * Sort the match collection by score type.
 */
void MatchCollection::sort(
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  )
{
  carp(CARP_DEBUG, "Sorting match collection.");

  // check if we are allowed to alter match_collection
  if(iterator_lock_){
    carp(CARP_FATAL,
         "Cannot sort a match collection when a match iterator is already"
         " instantiated");
  }

  // Switch to the equivalent sort key.
  SCORER_TYPE_T sort_by = NUMBER_SCORER_TYPES; // Initialize to nonsense.
  int (*compare_match_function)(const void*, const void*) 
    = (QSORT_COMPARE_METHOD)compareSp;
  switch (score_type) {
  case SP: 
    carp(CARP_DEBUG, "Sorting match collection by Sp.");
    sort_by = SP;
    compare_match_function = (QSORT_COMPARE_METHOD)compareSp;
    break;
  case EVALUE:
    sort_by = EVALUE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareEValue;
    break;
  case XCORR:
  case DECOY_XCORR_QVALUE:
  case LOGP_WEIBULL_XCORR: 
  case DECOY_XCORR_PEPTIDE_QVALUE:
  case DECOY_XCORR_PEP:
    carp(CARP_DEBUG, "Sorting match collection by XCorr.");
    sort_by = XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)compareXcorr;
    break;
  case TIDE_SEARCH_EXACT_PVAL:
    carp(CARP_DEBUG, "Sorting match collection by exact p-value.");
    sort_by = TIDE_SEARCH_EXACT_PVAL;
    compare_match_function = (QSORT_COMPARE_METHOD)compareExactPValue;
    break;  
  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR:
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
  case LOGP_WEIBULL_PEP:
    carp(CARP_DEBUG, "Sorting match collection by p-value.");
    sort_by = LOGP_BONF_WEIBULL_XCORR;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePValue;
    break;
/*
  case PERCOLATOR_SCORE:
    carp(CARP_INFO, "Sorting match collection by Percolator score.");
    sort_by = PERCOLATOR_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePercolatorScore;
    break;
*/
  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
  case PERCOLATOR_PEP:
    carp(CARP_DEBUG, "Sorting match collection by Percolator q-value.");
    sort_by = PERCOLATOR_QVALUE;
    compare_match_function = (QSORT_COMPARE_METHOD)comparePercolatorQValue;
    break;
/*
  case PERCOLATOR_SCAN:
    carp(CARP_INFO, "Sorting match collection by Percolator scan.");
    sort_by = PERCOLATOR_SCAN;
    compare_match_function = (QSORT_COMPARE_METHOD)compareSpectrum;
    break;
*/
  case QRANKER_SCORE:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker score.");
    sort_by = QRANKER_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareQRankerScore;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
  case QRANKER_PEP:
    carp(CARP_DEBUG, "Sorting match collection by Q-ranker q-value.");
    compare_match_function = (QSORT_COMPARE_METHOD)compareQRankerQValue;
    sort_by = QRANKER_QVALUE;
    break;

  case BARISTA_SCORE:
    carp(CARP_DEBUG, "Sorting match collection by barista score.");
    sort_by = BARISTA_SCORE;
    compare_match_function = (QSORT_COMPARE_METHOD)compareBaristaScore;
    break;

  case BARISTA_QVALUE:
  case BARISTA_PEPTIDE_QVALUE:
  case BARISTA_PEP:
    carp(CARP_DEBUG, "Sorting match collection by barista q-value.");
    compare_match_function = (QSORT_COMPARE_METHOD)compareBaristaQValue;
    sort_by = BARISTA_QVALUE;
    break;

  // Should never reach this point.
  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Something is terribly wrong in the sorting code!");
  }

  // Don't sort if it's already sorted.
  if (last_sorted_ == sort_by) {
    return;
  }

  // Do the sort.
  qsortMatch(match_, compare_match_function);
  last_sorted_ = sort_by;
}

/**
 * \brief Sort a match_collection by the given score type, grouping
 * matches by spectrum (if multiple spectra are present).
 */
void MatchCollection::spectrumSort(
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  ){

  // check if we are allowed to alter match_collection
  if(iterator_lock_){
    carp(CARP_FATAL,
         "Cannot alter match_collection when a match iterator is already"
         " instantiated");
  }

  switch(score_type){
  case SP: 
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumSp);
    last_sorted_ = SP;
    break;
  case TIDE_SEARCH_EXACT_PVAL:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareExactPValue);
    last_sorted_ = TIDE_SEARCH_EXACT_PVAL;    
    break;
  case XCORR:
  case LOGP_WEIBULL_XCORR: 
  case LOGP_BONF_WEIBULL_XCORR: 
  case LOGP_QVALUE_WEIBULL_XCORR: 
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
  case LOGP_WEIBULL_PEP:
  case DECOY_XCORR_QVALUE:
  case DECOY_XCORR_PEPTIDE_QVALUE:
  case DECOY_XCORR_PEP:
    /* If we are sorting on a per-spectrum basis, then the xcorr is
       good enough, even in the presence of p-values. */
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumXcorr);
    last_sorted_ = XCORR;
    break;

  case PERCOLATOR_SCORE:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumPercolatorScore);
    last_sorted_ = PERCOLATOR_SCORE;
    break;

  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
  case PERCOLATOR_PEP:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumPercolatorQValue);
    last_sorted_ = PERCOLATOR_QVALUE;
    break;
/*
  case PERCOLATOR_SCAN: 
    qsortMatch(
      match_,
      match_total_,
      (QSORT_COMPARE_METHOD)compareSpectrumScan
    );
    last_sorted_=PERCOLATOR_SCAN;
    break; 
*/
  case QRANKER_SCORE:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumQRankerScore);
    last_sorted_ = QRANKER_SCORE;
    break;

  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
  case QRANKER_PEP:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumQRankerQValue);
    last_sorted_ = QRANKER_QVALUE;
    break;

  case BARISTA_SCORE:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumBaristaScore);
    last_sorted_ = BARISTA_SCORE;
    break;

  case BARISTA_QVALUE:
  case BARISTA_PEPTIDE_QVALUE:
  case BARISTA_PEP:
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareSpectrumBaristaQValue);
    last_sorted_ = BARISTA_QVALUE;
    break;


  // Should never reach this point.
  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Something is terribly wrong in the sorting code!");
 }
}

/**
 * Assigns a rank for the given score type to each match.  First sorts
 * by the score type (if not already sorted).  Overwrites any existing
 * rank values, so it can be performed on a collection with matches
 * newly added to previously ranked matches.  Rank 1 is highest
 * score.  Matches with the same score will be given the same rank.
 *
 * \returns true, if populates the match rank in the match collection
 */
bool MatchCollection::populateMatchRank(
 SCORER_TYPE_T score_type ///< score type (SP, XCORR) by which to rank -in
 )
{
  carp(CARP_DETAILED_DEBUG, "Ranking matches by %i.", score_type);
  carp(CARP_DETAILED_DEBUG, "Collection currently ranked by %d", last_sorted_);
  // check if the match collection is in the correct sorted order
  if(last_sorted_ != score_type){
    // sort match collection by score type
    carp(CARP_DETAILED_DEBUG, "Sorting by score_type %i", score_type);
    sort(score_type);
  }

  // set match rank for all match objects that have been scored for
  // this type
  int cur_rank = 0;
  FLOAT_T cur_score = NOT_SCORED;
  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    FLOAT_T this_score = (*i)->getScore(score_type);
    
    if( NOT_SCORED == (*i)->getScore(score_type) ){
      char* seq = (*i)->getModSequenceStrWithMasses(MOD_MASS_ONLY);
      carp(CARP_WARNING, 
           "PSM spectrum %i charge %i sequence %s was NOT scored for type %i",
           (*i)->getSpectrum()->getFirstScan(),
           (*i)->getCharge(), seq,
           (int)score_type);
      free(seq);
    }

    // does this match have a higher score?
    if( this_score != cur_score ){
      cur_score = this_score;
      cur_rank++;
    }

    //    set_match_rank( cur_match, score_type, match_index+1);
    (*i)->setRank(score_type, cur_rank);

    carp(CARP_DETAILED_DEBUG, "Match rank %i, score %f", cur_rank, cur_score);
  }
  
  return true;
}

/**
 * For the #top_count ranked peptides, calculate the Weibull parameters
 *\returns true, if successfully calculates the Weibull parameters
 */
static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.0;       // For now, turn off the threshold.
static const FLOAT_T XCORR_SHIFT = 0.05;
static const FLOAT_T MIN_SP_SHIFT = -100.0;
static const FLOAT_T MAX_SP_SHIFT = 300.0;
static const FLOAT_T SP_SHIFT = 5.0;

/**
 * match_collection get, set method
 */

/**
 *\returns true, if the match collection has been scored by score_type
 */
bool MatchCollection::getScoredType(
  SCORER_TYPE_T score_type ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  )
{
  return scored_type_[score_type];
}

/**
 * sets the score_type to value
 */
void MatchCollection::setScoredType(
  SCORER_TYPE_T score_type, ///< the score_type (MATCH_SP, MATCH_XCORR) -in
  bool value
  )
{
  scored_type_[score_type] = value;
}

/**
 *
 */
void MatchCollection::getCustomScoreNames(
  vector<string>& custom_score_names
  ) {
  custom_score_names.clear();

  if (!match_.empty()) {
    match_.front()->getCustomScoreNames(custom_score_names);
  }

}

/**                                                                                                    
 * Set the filepath for all matches in the collection                                                  
 *\returns the associated file idx                                                                    
 */
int MatchCollection::setFilePath(
  const string& file_path  ///< File path to set                                                  
  ) {

  if (!match_.empty()) {
    vector<Crux::Match*>::iterator i = match_.begin();
    int file_idx = (*i)->setFilePath(file_path);
    for (i = i + 1; i != match_.end(); i++) {
      (*i)->setFileIndex(file_idx);
    }
    return file_idx;
  } else {
    carp(CARP_WARNING, "MatchCollection::setFilePath(): No matches in %s",file_path.c_str());
    return -1;
  }
}

/**
 *\returns true, if there is a  match_iterators instantiated by match collection 
 */
bool MatchCollection::getIteratorLock()
{
  return iterator_lock_;
}

/**
 *\returns the total match objects avaliable in current match_collection
 */
int MatchCollection::getMatchTotal()
{
  return match_.size();
}

bool MatchCollection::getHasDistinctMatches() {
  return has_distinct_matches_;
}

void MatchCollection::setHasDistinctMatches(bool distinct) {
  has_distinct_matches_ = distinct;
}


void MatchCollection::setExperimentSize(int size)
{
  experiment_size_ = size;
}

/**
 * \returns The total number of peptides searched for this spectrum,
 * target peptides for a target collection or decoy peptides for a
 * decoy collection.
 */
int MatchCollection::getExperimentSize()
{
  return experiment_size_;
}

/**
 * Sets the total number of target peptides searched for this
 * spectrum.  Only needs to be used by decoy match collections.
 */
void MatchCollection::setTargetExperimentSize(int numMatches){
  target_experiment_size_ = numMatches;
}

/**
 * \returns The number of target matches that this spectrum had.
 * Different than getExperimentSize() for decoy match collections.
 */
int MatchCollection::getTargetExperimentSize(){
  if( null_peptide_collection_ ){
    return target_experiment_size_;
  }
  return experiment_size_;
}

/**
 *\returns the charge of the spectrum that the match collection was created
 */
int MatchCollection::getCharge()
{
  return zstate_.getCharge();
}

/**
 * Must have been scored by Xcorr, returns error if not scored by Xcorr
 *\returns the delta cn value(difference in top and second ranked Xcorr values)
 */
FLOAT_T MatchCollection::getDeltaCn()
{
  // Check if xcorr value has been scored, thus delta cn value is valid
  if(scored_type_[XCORR]){
    return delta_cn_;
  }
  else{
    carp(CARP_ERROR, "must score match_collection with XCORR to get delta cn value");
    return 0.0;
  }
}

/**
 * \brief Transfer the Weibull distribution parameters, including the
 * correlation from one match_collection to another.  No check to see
 * that the parameters have been estimated.
 */
void MatchCollection::transferWeibull(
  MatchCollection* from_collection,
  MatchCollection* to_collection
  ){
  to_collection->eta_ = from_collection->eta_;
  to_collection->beta_ = from_collection->beta_;
  to_collection->shift_ = from_collection->shift_;
  to_collection->correlation_ = from_collection->correlation_;
}

/**
 * \brief Prints out the pepxml header to the output stream
 * passed in as a parameter.
 */

void MatchCollection::printXmlHeader(
  FILE* output
  ){
  if (output == NULL ){
    return;
  }
  time_t hold_time;
  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  char* enz_str = enzyme_type_to_string(enzyme);
  string database = get_string_parameter("protein-database");
  string msms_file = get_string_parameter("ms2 file");
  char* absolute_msms_path;
  if (msms_file.empty()) {
    absolute_msms_path = (char*) malloc(sizeof(char)*3);
    strcpy(absolute_msms_path, "NA");
  } else {
#if DARWIN
    char path_buffer[PATH_MAX];
    absolute_msms_path =  realpath(msms_file.c_str(), path_buffer);
#else
    absolute_msms_path =  realpath(msms_file.c_str(), NULL);
#endif
  }
  // Removes the extension from ms2 file path
  char* extension = strstr(absolute_msms_path, ".ms2");
  if (extension != NULL) (*extension) = '\0';
  
  
  MASS_TYPE_T isotopic_mass_type = get_mass_type_parameter("isotopic-mass");
  MASS_TYPE_T fragment_mass_type = get_mass_type_parameter("fragment-mass");

  const char* isotopic_mass;
  const char* fragment_mass;
  DIGEST_T digest = get_digest_type_parameter("digestion");
  int max_num_internal_cleavages;
  int min_number_termini;
  int missed_cleavage = get_int_parameter("missed-cleavages");
  if (missed_cleavage){
    max_num_internal_cleavages = get_int_parameter("max-length");
  } else {
    max_num_internal_cleavages = 0;
  }

  if (digest == FULL_DIGEST){
    min_number_termini = 2;
  } else if (digest == PARTIAL_DIGEST){
    min_number_termini = 1;
  } else {
    min_number_termini = 0;
  }
  
  if (isotopic_mass_type == AVERAGE){
    isotopic_mass = "average";
  } else {
    isotopic_mass = "monoisotopic";
  }

  if (fragment_mass_type == AVERAGE){
    fragment_mass =  "average";
  } else {
    fragment_mass =  "monoisotopic";
  }

  char* absolute_database_path = NULL;
  if (!database.empty()) {
#if DARWIN
    char path_buffer[PATH_MAX];
    absolute_database_path =  realpath(database.c_str(), path_buffer);
#else
    absolute_database_path =  realpath(database.c_str(), NULL);
#endif
  }
  
  hold_time = time(0);

  const char* xsi = "http://www.w3.org/2001/XMLSchema-instance";
  const char* xmlns = "http://regis-web.systemsbiology.net/pepXML";
  const char* schema_location = "/usr/local/tpp/schema/pepXML_v110.xsd";
  fprintf(output, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(output, "<?xml-stylesheet type=\"text/xsl\" href=\"\"?>\n");
  fprintf(output, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\""
          " xmlns:xsi=\"%s\" xsi:schemaLocation=\"%s %s\""
          " summary_xml=\"\">\n",
          ctime(&hold_time),
          xmlns, xsi, xmlns, schema_location);

  fprintf(output, "<msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" "
          "msModel=\"%s\" msIonization=\"%s\" msAnalyzer=\"%s\" "
          "msDectector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\" >\n",
          absolute_msms_path,
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA", // TODO, dummy value
          "NA" // TODO, dummy value
          );
  

  fprintf(output, "<sample_enzyme name=\"%s\">\n</sample_enzyme>\n", enz_str);

  fprintf(output, "<search_summary base_name=\"%s\" search_engine=\"%s\" "
          "precursor_mass_type=\"%s\" fragment_mass_type=\"%s\" "
          "out_data_type=\"%s\" out_data=\"%s\" search_id=\"%i\" >\n",
          absolute_msms_path,
          "Crux",
          isotopic_mass, // isotopic mass type is precursor mass type?
          fragment_mass,
          "NA", // TODO, dummy value
          "NA",
          1 // TODO, dummy value
          );

  
  fprintf(output, "<search_database local_path=\"%s\" type=\"%s\" />\n", 
          absolute_database_path, 
          "AA"
          );
  fprintf(output, "<enzymatic_search_constraint enzyme=\"%s\" "
          "max_num_internal_cleavages=\"%i\" min_number_termini=\"%i\"/>\n",
          enz_str,
          max_num_internal_cleavages,
          min_number_termini
          );

#ifndef DARWIN
  free(absolute_msms_path);
  free(absolute_database_path);
#endif
  free(enz_str);


  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A'+ ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");
  int aa = 0;

  // static amino acid modifications
  for (aa = (int)'A'; aa < alphabet_size-1; aa++){
    aa_str[0] = (char)aa;
    double mod = get_double_parameter(aa_str);
    double mass = get_mass_amino_acid(aa, isotopic_type);
    
    if (mod != 0 ){
      fprintf(output, "<aminoacid_modification aminoacid=\"%s\" mass=\"%f\" "
              "massdiff=\"%f\" variable=\"%s\" />\n",
              aa_str,
              mass + mod,
              mod,
              "N" // N if static modification
              );      
    }
  }
  
  // variable amino acid modifications
  AA_MOD_T** mod_list = NULL;
  int num_mods = get_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    
    bool* aas_modified = aa_mod_get_aa_list(mod_list[mod_idx]);
    for (int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
      if (aas_modified[aa_idx] == true ){
        int aa = (aa_idx+'A');
        FLOAT_T aa_mass = get_mass_amino_acid(aa , isotopic_type);
        fprintf(output, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" "
                "massdiff=\"%f\" variable=\"%s\" />\n",
                aa,
                aa_mass + mod_mass,
                mod_mass,
                "Y" // Y if variable modification
                );    
      }
    }

  }

  // terminal modifciations
  // variable
  num_mods = get_c_mod_list(&mod_list); // variable c mods
  for(int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    fprintf(output, "<terminal_modification terminus=\"c\" "
            "mass=\"%f\" massdiff=\"%f\" variable=\"Y\" />\n",
            MASS_OH + mod_mass,
            mod_mass
            );
  }
  num_mods = get_n_mod_list(&mod_list); // variable n mods
  for(int mod_idx = 0; mod_idx < num_mods; mod_idx++){
    FLOAT_T mod_mass = aa_mod_get_mass_change(mod_list[mod_idx]);
    fprintf(output, "<terminal_modification terminus=\"n\" "
            "mass=\"%f\" massdiff=\"%f\" variable=\"Y\" />\n",
            MASS_H_MONO + mod_mass,
            mod_mass
            );
  }
  // fixed
  if( get_num_fixed_mods() != 0 ){
    get_all_aa_mod_list(&mod_list);
    int fixed_mod_idx = get_fixed_mod_index(N_TERM); // fixed n mods
    if( fixed_mod_idx > -1 ){
      fprintf(output, "<terminal_modification terminus=\"n\" "
              "mass=\"?\" massdiff=\"%f\" variable=\"N\" />\n",
              aa_mod_get_mass_change(mod_list[fixed_mod_idx])
              );
    }

    fixed_mod_idx = get_fixed_mod_index(C_TERM); // fixed c mods
    if( fixed_mod_idx > -1 ){
      fprintf(output, "<terminal_modification terminus=\"c\" "
              "mass=\"?\" massdiff=\"%f\" variable=\"N\" />\n",
              aa_mod_get_mass_change(mod_list[fixed_mod_idx])
              );
    }

  }


  print_parameters_xml(output);
  
  fprintf(output, "</search_summary>\n");

}

/**
 * Write header for .sqt file.  Assumes only sequest-search is writing
 * this file type.
 */
void MatchCollection::printSqtHeader(
 FILE* output, 
 const char* type, 
 string database,
 int num_proteins,
 bool exact_pval_search){  
  if( output == NULL ){
    return;
  }

  time_t hold_time;
  hold_time = time(0);

  bool decoy = false;
  if( strcmp(type, "decoy") == 0 ){
    decoy = true;
  }

  fprintf(output, "H\tSQTGenerator Crux\n");
  fprintf(output, "H\tSQTGeneratorVersion 1.0\n");
  fprintf(output, "H\tComment Crux was written by...\n");
  fprintf(output, "H\tComment ref...\n");
  fprintf(output, "H\tStartTime\t%s", ctime(&hold_time));
  fprintf(output, "H\tEndTime                               \n");

  fprintf(output, "H\tDatabase\t%s\n", database.c_str());

  if(decoy){
  fprintf(output, "H\tComment\tDatabase shuffled; these are decoy matches\n");
  }
  fprintf(output, "H\tDBSeqLength\t?\n");
  fprintf(output, "H\tDBLocusCount\t%d\n", num_proteins);

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");
  char temp_str[64];
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tPrecursorMasses\t%s\n", temp_str);
  
  mass_type = get_mass_type_parameter("fragment-mass");
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tFragmentMasses\t%s\n", temp_str); //?????????

  double tol = get_double_parameter("precursor-window");
  fprintf(output, "H\tAlg-PreMasTol\t%.1f\n",tol);
  fprintf(output, "H\tAlg-FragMassTol\t%.2f\n", 
          get_double_parameter("mz-bin-width") / 2.0);
  fprintf(output, "H\tAlg-XCorrMode\t0\n");

  SCORER_TYPE_T score = get_scorer_type_parameter("prelim-score-type");
  fprintf(output, "H\tComment\tpreliminary algorithm %s\n", 
          scorer_type_to_string(score));

  score = get_scorer_type_parameter("score-type");
  fprintf(output, "H\tComment\tfinal algorithm %s\n", 
          scorer_type_to_string(score));

  int aa = 0;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = get_mass_type_parameter("isotopic-mass");

  for(aa = (int)'A'; aa < alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    double mod = get_double_parameter(aa_str);
    if( mod != 0 ){
      //      double mass = mod + get_mass_amino_acid(aa, isotopic_type);
      double mass = get_mass_amino_acid(aa, isotopic_type);
      fprintf(output, "H\tStaticMod\t%s=%.3f\n", aa_str, mass);
    }
  }
  // print dynamic mods, if any
  // format DiffMod <AAs><symbol>=<mass change>
  AA_MOD_T** aa_mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char* aa_list_str = aa_mod_get_aa_list_string(aamod);
    char aa_symbol = aa_mod_get_symbol(aamod);
    double mass_dif = aa_mod_get_mass_change(aamod);

    fprintf(output, "H\tDiffMod\t%s%c=%+.2f\n", aa_list_str, 
            aa_symbol, mass_dif);
    free(aa_list_str);
  }
  num_mods = get_c_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    fprintf(output, "H\tComment\tMod %c is a C-terminal modification\n",
            aa_symbol);
  }

  num_mods = get_n_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aa_mod_get_symbol(aamod);

    fprintf(output, "H\tComment\tMod %c is a N-terminal modification\n",
            aa_symbol);
  }



  //for letters in alphabet
  //  double mod = get_double_parameter(letter);
  //  if mod != 0
  //     double mass = mod + getmass(letter);
  //     fprintf(output, "H\tStaticMod\t%s=%.3f\n", letter, mass);
  //  fprintf(output, "H\tStaticMod\tC=160.139\n");
  fprintf(output, "H\tAlg-DisplayTop\t%d\n", 
          //          get_int_parameter("max-sqt-result")); 
          get_int_parameter("top-match")); 
  // this is not correct for an sqt from analzyed matches

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  char* enz_str = enzyme_type_to_string(enzyme);
  char* dig_str = digest_type_to_string(digestion);
  char custom_str[SMALL_BUFFER];
  if( enzyme == CUSTOM_ENZYME){
    string rule = get_string_parameter("custom-enzyme");
    sprintf(custom_str, ", custom pattern: %s", rule.c_str());
  }else{
    custom_str[0] = 0;
  }
  fprintf(output, "H\tEnzymeSpec\t%s-%s%s\n", enz_str, dig_str, custom_str);
  free(enz_str);
  free(dig_str);
  // write a comment that says what the scores are
  fprintf(output, "H\tLine fields: S, scan number, scan number, "
          "charge, 0, server, experimental mass, total ion intensity, "
          "lowest Sp, number of matches\n");
  if (exact_pval_search) {
    fprintf(output, "H\tLine fields: M, rank by xcorr score, rank by sp score, "
            "peptide mass, deltaCn, exact P-value, recalibrated xcorr, sp score, number ions matched, "
            "total ions compared, sequence, validation status\n");
  } else {
    fprintf(output, "H\tLine fields: M, rank by xcorr score, rank by sp score, "
            "peptide mass, deltaCn, xcorr score, sp score, number ions matched, "
            "total ions compared, sequence, validation status\n");
  }
}

/**
 * Print the header line for a tab-delimited file.
 */
void MatchCollection::printTabHeader(FILE* output){

  if( output == NULL ){
    return;
  }

  int column_idx;
  for (column_idx = 0; column_idx < NUMBER_MATCH_COLUMNS; column_idx++) {
    fprintf(output, "%s", get_column_header(column_idx));
    if (column_idx < NUMBER_MATCH_COLUMNS - 1) {
      fprintf(output, "\t");
    } else {
      fprintf(output, "\n");
    }
  }
}


/**
 * Print Footer lines for xml files
 */
void MatchCollection::printXmlFooter(FILE* output){
  if (output == NULL ){
    return;
  }
  
  fprintf(output, "</msms_run_summary>\n");
  fprintf(output, "</msms_pipeline_analysis>\n");
}


/**
 * \brief Print the given match collection for several spectra to
 * xml files only. Takes the spectrum information from the
 * matches in the collection. At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter. 
 */
void MatchCollection::printMultiSpectraXml(
  PepXMLWriter* output
){
  carp(CARP_DETAILED_DEBUG, "Writing matches to xml file");
  // reuse these for each match
  vector<string> protein_ids;
  vector<string> protein_descriptions;
  double* scores = new double[NUMBER_SCORER_TYPES];

  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    Match* cur_match = *i;
    bool is_decoy = cur_match->getNullPeptide();
    Spectrum* spectrum = cur_match->getSpectrum();
    double cur_ln_experiment_size=0;
    if (! is_decoy){
      int* ranks =new int[NUMBER_SCORER_TYPES];
      ranks[XCORR]=-1; 
      if( scored_type_[XCORR] ){
        ranks[XCORR] = cur_match->getRank(XCORR);
      }else if(scored_type_[SP]){
        ranks[SP]=cur_match->getRank(SP);
        scores[SP]=cur_match->getScore(SP);
      }
      char* peptide_sequence = cur_match->getSequence();
      char* mod_peptide_sequence = cur_match->getModSequenceStrWithMasses(
                             get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = cur_match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, 
                                                 protein_descriptions);
      for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
        if( scored_type_[score_idx] == true ){
          scores[score_idx] = cur_match->getScore((SCORER_TYPE_T)score_idx);
          ranks[score_idx]=cur_match->getRank((SCORER_TYPE_T)score_idx);
         
        }
  
      }
      unsigned num_matches = getTargetExperimentSize(); 
      if(isDecoy())
        num_matches=getExperimentSize();
      else 
        num_matches=getTargetExperimentSize(); 
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        cur_match->getNeutralMass(),
        cur_match->getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->getPeptideMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        cur_match->getDeltaCn(),
        scored_type_,
        scores,
        cur_match->getBYIonMatched(),
        cur_match->getBYIonPossible(),
        num_matches
      );
    }
  }
  
}

/**
 * \brief Print the psm features to file in xml format
 *
 * Prints a spectrum_query tag which encompasses the search_hit tag
 * which represents peptide to spectra match.
 *
 * returns true, if succesfully printed xml format of PSMs, else false
 *
 */

bool MatchCollection::printXml(
  PepXMLWriter* output,
  int top_match,
  Spectrum* spectrum,
  SCORER_TYPE_T main_score
  )
{
  if ( output == NULL || spectrum == NULL || match_.empty()){
    return false;
  }

  // calculate delta_cn and populate fields in the matches
  calculateDeltaCn();

  // for deciding when to quit
  int count = 0;
  int last_rank = 0;

  // reuse these for each match
  vector<string> protein_ids;
  vector<string> protein_descriptions;
  bool* scores_computed = new bool[NUMBER_SCORER_TYPES];
  
  for(int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
    scores_computed[score_idx] = false;
  }
  scores_computed[main_score] = true;
  if( scored_type_[SP]) {
    scores_computed[SP] = true;
  }
  scores_computed[TIDE_SEARCH_EXACT_PVAL] = exact_pval_search_;
  scores_computed[TIDE_SEARCH_REFACTORED_XCORR] = exact_pval_search_;
  scored_type_[TIDE_SEARCH_EXACT_PVAL] = exact_pval_search_;
  scored_type_[TIDE_SEARCH_REFACTORED_XCORR] = exact_pval_search_;
  scores_computed[main_score] = !exact_pval_search_;

  double* scores = new double[NUMBER_SCORER_TYPES];
  int* ranks=new int[NUMBER_SCORER_TYPES];

  Match* match = NULL;
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator = 
    new MatchIterator(this, main_score, true);
  // iterate over matches
  while(match_iterator->hasNext()){
    match = match_iterator->next();
    int cur_rank = match->getRank(main_score);   
    if(scored_type_[XCORR])
      ranks[XCORR]=match->getRank(XCORR);
    if(scored_type_[SP]){
      ranks[SP]= match->getRank(SP);
      scores[SP]= match->getScore(SP);
    }
    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){
      
      char* peptide_sequence = match->getSequence();
      char* mod_peptide_sequence = match->getModSequenceStrWithMasses(
                            get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, 
                                                 protein_descriptions);
     for(int score_idx=0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
      if(scored_type_[score_idx])
        scores[score_idx] = match->getScore((SCORER_TYPE_T)score_idx);
      
     }   
     unsigned num_matches= getTargetExperimentSize(); 
     if(isDecoy())
       num_matches= getExperimentSize(); 
     else
       num_matches= getTargetExperimentSize(); 
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        zstate_.getNeutralMass(),
        zstate_.getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->getPeptideMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        match->getDeltaCn(),
        scores_computed,
        scores,
        match->getBYIonMatched(),
        match->getBYIonPossible(), 
        num_matches
      );
      count++;
      last_rank = cur_rank;
      free(peptide_sequence);
      free(mod_peptide_sequence);
      free(flanking_aas);
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d xml matches", 
       count, match_.size());

  delete match_iterator;
  delete scores_computed;
  delete scores;
  delete ranks;

  return true;
}


/**
 * \brief Print the psm features to file in sqt format.
 *
 * Prints one S line, 'top_match' M lines, and one locus line for each
 * peptide source of each M line.
 * Assumes one spectrum per match collection.  Only crux
 * sequset-search produces sqt files so the two scores are always Sp
 * and xcorr.
 * Possible side effects: Calculates delta_cn and stores in each match
 *\returns true, if sucessfully print sqt format of the PSMs, else false 
 */
bool MatchCollection::printSqt(
  FILE* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  Spectrum* spectrum           ///< the spectrum to print sqt -in
  )
{

  if( output == NULL || spectrum == NULL || match_.empty() ){
    return false;
  }

  SpectrumZState& zstate = zstate_; 
  int num_matches = experiment_size_;

  // calculate delta_cn and populate fields in the matches
  calculateDeltaCn();

  // First, print spectrum info
  spectrum->printSqt(output, num_matches, zstate);
  
  Match* match = NULL;
  
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator =
    new MatchIterator(this, XCORR, true);
  
  // Second, iterate over matches, prints M and L lines
  while(match_iterator->hasNext()){
    match = match_iterator->next();    

    // print only up to max_rank_result of the matches
    if( match->getRank(XCORR) > top_match ){
      break;
    }// else

    match->printSqt(output);

  }// next match
  
  // make sure top_scoring_sp_ has been set
  if( top_scoring_sp_ == NULL){
    carp(CARP_DEBUG, "Top scoring SP was not set.");
  } else if( top_scoring_sp_->getRank(XCORR) > top_match ){
    // print the match with Sp rank==1 if its xcorr rank > top_match rank.  
    top_scoring_sp_->printSqt(output);
  }
  
  delete match_iterator;
  
  return true;
}

/**
 * \brief Print the psm features to file in tab delimited format.
 *
 * Matches will be sorted by main_score and the ranks of those scores
 * will be used to determine how many matches are printed for each
 * spectrum.
 * \returns true, if sucessfully print tab-delimited format of the
 * PSMs, else false  
 */
bool MatchCollection::printTabDelimited(
  MatchFileWriter* output,                  ///< the output file -out
  int top_match,                 ///< the top matches to output -in
  Spectrum* spectrum,          ///< the spectrum to print sqt -in
  SCORER_TYPE_T main_score       ///< the main score to report -in
  )
{

  if( output == NULL || spectrum == NULL || match_.empty() ){
    return false;
  }
  int num_target_matches = getTargetExperimentSize();
  int num_decoy_matches = getExperimentSize();
  int scan_num = spectrum->getFirstScan();
  FLOAT_T spectrum_precursor_mz = spectrum->getPrecursorMz();

  Match* match = NULL;
  
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator = 
    new MatchIterator(this, main_score, true);
  int count = 0;
  int last_rank = 0;

  // iterate over matches
  while(match_iterator->hasNext()){
    match = match_iterator->next();
    int cur_rank = match->getRank(main_score);

    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if( count < top_match || last_rank == cur_rank ){

      match->printTab(this, output, scan_num, 
                      spectrum_precursor_mz, 
                      num_target_matches, num_decoy_matches);
      count++;
      last_rank = cur_rank;
    } else if( count >= top_match && last_rank != cur_rank ) {
      break;
    } // else see if there is one more tie to print

  }// next match
  
  carp(CARP_DETAILED_DEBUG, "printed %d out of %d tab matches", 
       count, num_target_matches);

  delete match_iterator;
  
  return true;
}

/**
 * Retrieve the calibration parameter eta.
 */
FLOAT_T MatchCollection::getCalibrationEta()
{
  return eta_;
}

/**
 * Retrieve the calibration parameter beta.
 */
FLOAT_T MatchCollection::getCalibrationBeta()
{
  return beta_;
}

/**
 * Retrieve the calibration parameter shift.
 */
FLOAT_T MatchCollection::getCalibrationShift()
{
  return shift_;
}

/**
 * Retrieve the calibration correlation.
 */
FLOAT_T MatchCollection::getCalibrationCorr()
{
  return correlation_;
}

/**
 * \brief Print the given match collection for several spectra to
 * tab-delimited files only.  Takes the spectrum information from the
 * matches in the collection.  At least for now, prints all matches in
 * the collection rather than limiting by top-match parameter.  Uses
 * SP as preliminary score and XCORR as main score.
 */
void MatchCollection::printMultiSpectra(
 MatchFileWriter* tab_file, 
 MatchFileWriter* decoy_tab_file
  ){

  carp(CARP_DETAILED_DEBUG, "Writing matches to file");

  // if file location is target (i.e. tdc=T), print all to target
  MatchFileWriter* decoy_file = decoy_tab_file;
  if( get_boolean_parameter("tdc") == true ){
    decoy_file = tab_file;
  }

  // for each match, get spectrum info, determine if decoy, print
  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    Match* cur_match = *i;
    bool is_decoy = cur_match->getNullPeptide();
    Spectrum* spectrum = cur_match->getSpectrum();
    int scan_num = spectrum->getFirstScan();
    FLOAT_T mz = spectrum->getPrecursorMz();
    FLOAT_T num_psm_per_spec = cur_match->getLnExperimentSize();
    num_psm_per_spec = expf(num_psm_per_spec) + 0.5; // round to nearest int
    int num_target_psm_per_spec = cur_match->getTargetExperimentSize();

    if( is_decoy ){
      cur_match->printTab(this, decoy_file, scan_num, mz, 
                          (int)num_psm_per_spec, num_target_psm_per_spec);
    }
    else{
      cur_match->printTab(this, tab_file, scan_num, mz,
                          (int)num_psm_per_spec, num_target_psm_per_spec);
    }

  }

}

/*******************************************
 * match_collection post_process extension
 ******************************************/

/**
 * Create a list of file names that we will look for in the psm directory.
 */
void set_possible_names(vector<string>& possible_names, SET_TYPE_T type){

  // if a specific file has been requested, return just that file name
  string psm_filename = get_string_parameter("input PSMs");
  if (!psm_filename.empty()){
    possible_names.push_back(psm_filename);
    return;
  }

  // else, decide on the search string for targets or decoys
  switch(type){
  case SET_TARGET:
    possible_names.push_back("search.target.txt");
    possible_names.push_back("sequest.target.txt");
    break;
  case SET_DECOY1:
    possible_names.push_back("search.decoy.txt");
    possible_names.push_back("search.decoy-1.txt");
    possible_names.push_back("sequest.decoy.txt");
    possible_names.push_back("sequest.decoy-1.txt");
    break;
  case SET_DECOY2:
    possible_names.push_back("search.decoy-2.txt");
    possible_names.push_back("sequest.decoy-2.txt");
    break;
  case SET_DECOY3:
    possible_names.push_back("search.decoy-3.txt");
    possible_names.push_back("sequest.decoy-3.txt");
    break;
  }
}

/**
 * Read files in the directory and return the names of target or
 * decoy files to use for post-search commands.
 * \returns Vector parameter filled with names of target or decoy
 * files.
 */
void get_target_decoy_filenames(vector<string>& target_decoy_names,
                                DIR* directory,
                                SET_TYPE_T type){
  if( directory == NULL ){
    carp(CARP_FATAL, "Cannot read files from NULL directory.");
  }

  // first see if there is a specific file to open
  string psm_file = get_string_parameter("input PSMs");
  if (!psm_file.empty()) {
    // strip off the path
    char** name_path = parse_filename_path(psm_file);
    target_decoy_names.push_back(name_path[0]);
    free(name_path[0]);
    free(name_path[1]);
    free(name_path);
    return;
  }

  // look for both files from search-for-matches and sequest-search
  vector<string> possible_names;

  set_possible_names(possible_names, type);

  // open the directory
  struct dirent* directory_entry = NULL;
  rewinddir(directory);
  // for each file, compare to each name, if it matches, add to the list
  while((directory_entry = readdir(directory))){
    for(int name_idx = 0; name_idx < (int)possible_names.size(); name_idx++){
      string filename = directory_entry->d_name;
      if( filename.find(possible_names[name_idx]) != string::npos ){
        target_decoy_names.push_back(filename);
      } 
    }
  }

  // check that files are only sequest or search, not both
  bool found_search = false;
  bool found_sequest = false;
  if( target_decoy_names.size() > 1 ){
    for(int name_idx = 0; name_idx < (int)target_decoy_names.size();name_idx++){
      // don't look for just "search" and "sequest" in case they are in
      // the fileroot
      if( target_decoy_names[name_idx].find(possible_names.front()) 
          != string::npos ){
        found_search = true;
      }
      if( target_decoy_names[name_idx].find(possible_names.back()) 
          != string::npos ){
        found_sequest = true;
      }
    }
  }

  if( found_search && found_sequest ){
    carp(CARP_FATAL, "Cannot analyze results from both crux search-for-matches "
         " and sequest-search.  Please remove one from the directory.");
  }
  // check that headers are all the same??
}


/**
 * parse all the match objects and add to match collection
 *\returns true, if successfully parse all PSMs in result_file, else false
 */
bool MatchCollection::extendTabDelimited(
  Database* database, ///< the database holding the peptides -in
  MatchFileReader& result_file,   ///< the result file to parse PSMs -in
  Database* decoy_database ///< the database holding the decoy peptides -in
  )
{
  Match* match = NULL;

  FLOAT_T delta_cn = 0;
  FLOAT_T ln_delta_cn = 0;
  FLOAT_T ln_experiment_size = 0;

  // only for post_process_collections
  if(!post_process_collection_){
    carp(CARP_ERROR, "Must be a post process match collection to extend.");
    return false;
  }

  while (result_file.hasNext()) {

    /*** get spectrum specific features ***/
    zstate_.setNeutralMass(
      result_file.getFloat(SPECTRUM_NEUTRAL_MASS_COL),
      result_file.getInteger(CHARGE_COL));
    scored_type_[DELTA_CN] = scored_type_[DELTA_CN] || !result_file.empty(DELTA_CN_COL);
    delta_cn = result_file.getFloat(DELTA_CN_COL);
    if (delta_cn <= 0.0) {
      ln_delta_cn = 0;
    } else {
      ln_delta_cn = logf(delta_cn);
    }
    if (!result_file.empty(DISTINCT_MATCHES_SPECTRUM_COL)) {
      has_distinct_matches_ = true;
      ln_experiment_size = log(result_file.getFloat(DISTINCT_MATCHES_SPECTRUM_COL));
    } else if (!result_file.empty(MATCHES_SPECTRUM_COL)) {
      ln_experiment_size = log(result_file.getFloat(MATCHES_SPECTRUM_COL));
    } else {
      ln_experiment_size = 0;
    }

    //TODO: Parse all boolean indicators for scores
    scored_type_[SP] = !result_file.empty(SP_SCORE_COL);

    scored_type_[XCORR] = !result_file.empty(XCORR_SCORE_COL);

    scored_type_[TIDE_SEARCH_EXACT_PVAL] = !result_file.empty(EXACT_PVALUE_COL);

    scored_type_[TIDE_SEARCH_REFACTORED_XCORR] = !result_file.empty(REFACTORED_SCORE_COL);
    
    scored_type_[EVALUE] = !result_file.empty(EVALUE_COL);
    
    scored_type_[DECOY_XCORR_QVALUE] = !result_file.empty(DECOY_XCORR_QVALUE_COL);

/* TODO
    match_collection -> 
      scored_type[LOGP_WEIBULL_XCORR] = 
      result_file.getString("logp weibull xcorr") != "";
*/

    scored_type_[LOGP_BONF_WEIBULL_XCORR] = 
      !result_file.empty(PVALUE_COL);

    scored_type_[PERCOLATOR_QVALUE] = 
      !result_file.empty(PERCOLATOR_QVALUE_COL);

    scored_type_[PERCOLATOR_SCORE] = 
      !result_file.empty(PERCOLATOR_SCORE_COL);

    scored_type_[LOGP_QVALUE_WEIBULL_XCORR] = 
      !result_file.empty(WEIBULL_QVALUE_COL);
  
    scored_type_[QRANKER_SCORE] = 
      !result_file.empty(QRANKER_SCORE_COL);
    
    scored_type_[QRANKER_QVALUE] = 
      !result_file.empty(QRANKER_QVALUE_COL);

    scored_type_[BARISTA_SCORE] =
      !result_file.empty(BARISTA_SCORE_COL);

    scored_type_[BARISTA_QVALUE] =
      !result_file.empty(BARISTA_QVALUE_COL);

    post_scored_type_set_ = true;

    // parse match object
    match = Match::parseTabDelimited(result_file, database, decoy_database);
    if (match == NULL) {
      carp(CARP_ERROR, "Failed to parse tab-delimited PSM match");
      return false;
    }

    //set all spectrum specific features to parsed match
    match->setZState(zstate_);
    match->setDeltaCn(delta_cn);
    match->setDeltaLCn(ln_delta_cn);
    match->setLnExperimentSize(ln_experiment_size);    
    //add match to match collection.
    addMatchToPostMatchCollection(match);
    //increment pointer.
    result_file.next();
  }

  return true;
}



/**
 * \brief Adds the match to match_collection by copying the pointer.
 * 
 * \returns true if successfully adds the match to the
 * match_collection, else false 
 */
bool MatchCollection::addMatch(
  Match* match ///< the match to add -in
  )
{
  if( match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match.");
  }

  // add a new match to array
  match_.push_back(match);
  match->incrementPointerCount();
  
  return true;
}

/**
 * Adds the match object to match_collection
 * Only for post_process (i.e. post search) match_collections.  Keeps
 * track of all peptides in a hash table.
 * \returns true if successfully adds the match to the
 * match_collection, else false 
 */
// this method renamed so that a more general add_match_to_match_collection could be implemented
bool MatchCollection::addMatchToPostMatchCollection(
  Match* match ///< the match to add -in
  )
{
  if( match == NULL ){
    carp(CARP_FATAL, "Cannot add NULL match to NULL collection.");
  }

  // only for post_process_collections
  if(!post_process_collection_){
    carp(CARP_ERROR, "Must be a post process match collection to add a match.");
    return false;
  }

  // add a new match to array
  match_.push_back(match);
  match->incrementPointerCount();
  
  // DEBUG, print total peptided scored so far
  if(match_.size() % 1000 == 0){
    carp(CARP_INFO, "parsed PSM: %d", match_.size());
  }
  
  return true;
}

/**
 * Fill the match objects score with the given the array, and populate
 * the corresponding ranks.  The match object order must not have been
 * altered since scoring.  The result array size must equal the number
 * of matches in the given match collection.  After the function
 * completes, the match collection is sorted by the specified score,
 * unless preserve_order is set to true.
 */
void MatchCollection::fillResult(
  double* results,           ///< array of scores -in
  SCORER_TYPE_T score_type,  ///< The score type of the results to fill -in
  bool preserve_order   ///< preserve match order?
  )
{
  Match** match_array = NULL;
  SCORER_TYPE_T score_type_old = last_sorted_;

  // iterate over match object in collection, set scores
  int match_idx = 0;
  for(; match_idx < match_.size(); ++match_idx){
    Match* match = match_[match_idx];
    match->setScore(score_type, results[match_idx]);    
  }
  
  // if need to preserve order store a copy of array in original order 
  if(preserve_order){
    match_array = (Match**)mycalloc(match_.size(), sizeof(Match*));
    for(match_idx=0; match_idx < match_.size(); ++match_idx){
      match_array[match_idx] = match_[match_idx];
    }
  }

  // populate the rank of match_collection
  if(!populateMatchRank(score_type)){
    carp(CARP_FATAL, "failed to populate match rank in match_collection");
  }
  
  // restore match order.
  if(preserve_order){
    for(match_idx=0; match_idx < match_.size(); ++match_idx){
      match_[match_idx] = match_array[match_idx];
    }
    last_sorted_ = score_type_old;
    free(match_array);
  }

  scored_type_[score_type] = true;
}

/**
 * \brief Calculate the delta_cn of each match and populate the field.
 * 
 * Delta_cn is the normalized difference between xcorrs of different ranks.
 * match[i] = (match[i] - match[i+1]) / match[i].
 * Sorts match_collection by xcorr, if necessary.
 * 
 */
bool MatchCollection::calculateDeltaCn(){

  if( scored_type_[XCORR] == false ){
    carp(CARP_WARNING, 
      "Delta_cn not calculated because match collection not scored for xcorr");
    return false;
  }

  // sort, if not already
  // N.B. Can't use sort_match_collection because iterator already exists!
  if( last_sorted_ != XCORR ){
    qsortMatch(match_, (QSORT_COMPARE_METHOD)compareXcorr);
    last_sorted_ = XCORR;
  }

  if (match_.size() < 2) {
    return true;
  }

  FLOAT_T last_xcorr=0.0;
  FLOAT_T delta_cn = 0.0;
  FLOAT_T delta_lcn = 0.0;
  FLOAT_T next_xcorr=0.0;
  FLOAT_T current_xcorr = 0 ; 
  last_xcorr = match_.back()->getScore(XCORR);
  for (size_t idx = 0 ;idx < match_.size();idx++) { 
    current_xcorr = match_[idx]->getScore(XCORR);
    if (idx+1<=match_.size()-1)
      next_xcorr=match_[idx+1]->getScore(XCORR);
    delta_cn = (current_xcorr - next_xcorr) / max(current_xcorr, (FLOAT_T)1.0);
    delta_lcn = (current_xcorr - last_xcorr) / max(current_xcorr, (FLOAT_T)1.0);
  
    if(fabs(delta_cn)== numeric_limits<FLOAT_T>::infinity()){
      carp(CARP_DEBUG, "delta_cn was %f and set to zero. XCorr score is %f", delta_cn, current_xcorr);
      delta_cn = 0.0;
    }   
    if(fabs(delta_lcn) == numeric_limits<FLOAT_T>::infinity()){
      carp(CARP_DEBUG, "delta_lcn was %f and set to zero. XCorr score is %f", delta_lcn, current_xcorr);
      delta_lcn = 0.0;
    }   
    match_[idx]->setDeltaCn(delta_cn);
    match_[idx]->setDeltaLCn(delta_lcn);
  }   

  return true;
}


/**********************************
 * match_collection get, set methods
 **********************************/

/**
 * \returns true if the match_collection only contains decoy matches,
 * else (all target or mixed) returns false.
 */
bool MatchCollection::isDecoy()
{
  return null_peptide_collection_;
}

/**
 * Try setting the match collection's charge.  Successful if the
 * current charge is 0 (i.e. hasn't yet been set) or if the current
 * charge is the same as the given value.  Otherwise, returns false
 *
 * \returns true if the match_collection's charge state was changed.
 */

bool MatchCollection::setZState(
  SpectrumZState& zstate ///< new zstate
  ) {

  if (getCharge() == 0) {
    zstate_ = zstate;
    return true;
  } else {
    //error
    carp(CARP_WARNING, "Cannot change the zstate of a match collection "
        "once it has been set.");
    return false;
  }
}

// cheater functions for testing

void MatchCollection::forceScoredBy(SCORER_TYPE_T type){
  scored_type_[type] = true;
}

/**
 * Extract a given type of score into an array.  The array is
 * allocated here and must be freed by the caller.
 */
FLOAT_T* MatchCollection::extractScores(
  SCORER_TYPE_T       score_type ///< Type of score to extract.
)
{
  FLOAT_T* return_value = (FLOAT_T*)mycalloc(match_.size(),
                                             sizeof(FLOAT_T));

  MatchIterator* match_iterator =
    new MatchIterator(this, score_type, false);
  int idx = 0;
  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    return_value[idx] = match->getScore(score_type);
    idx++;
  }
  delete match_iterator;

  return(return_value);
}

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
void MatchCollection::assignQValues(
  const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
  SCORER_TYPE_T score_type
){

  // Iterate over the matches filling in the q-values
  MatchIterator* match_iterator = 
    new MatchIterator(this, score_type, false);

  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    FLOAT_T score = match->getScore(score_type);

    // Retrieve the corresponding q-value.
    map<FLOAT_T, FLOAT_T>::const_iterator map_position 
      = score_to_qvalue_hash->find(score);
    if (map_position == score_to_qvalue_hash->end()) {
      carp(CARP_FATAL,
           "Cannot find q-value corresponding to score of %g.",
           score);
    }
    FLOAT_T qvalue = map_position->second;

    /* If we're given a base score, then store the q-value.  If we're
       given a q-value, then store the peptide-level q-value. */
    SCORER_TYPE_T derived_score_type = INVALID_SCORER_TYPE;
    switch (score_type) {
    case XCORR:
    case TIDE_SEARCH_EXACT_PVAL:
      derived_score_type = DECOY_XCORR_QVALUE;
      break;
    case DECOY_XCORR_QVALUE:
      derived_score_type = DECOY_XCORR_PEPTIDE_QVALUE;
      break;
    case EVALUE:
      derived_score_type = DECOY_EVALUE_QVALUE;
      break;
    case DECOY_EVALUE_QVALUE:
      derived_score_type = DECOY_EVALUE_PEPTIDE_QVALUE;
      break;
    case LOGP_BONF_WEIBULL_XCORR: 
      derived_score_type = LOGP_QVALUE_WEIBULL_XCORR;
      break;
    case LOGP_QVALUE_WEIBULL_XCORR:
      derived_score_type = LOGP_PEPTIDE_QVALUE_WEIBULL;
      break;
    case PERCOLATOR_SCORE:
      derived_score_type = PERCOLATOR_QVALUE;
      break;
    case PERCOLATOR_QVALUE:
      derived_score_type = PERCOLATOR_PEPTIDE_QVALUE;
      break;
    case QRANKER_SCORE:
      derived_score_type = QRANKER_QVALUE;
      break;
    case QRANKER_QVALUE:
      derived_score_type = QRANKER_PEPTIDE_QVALUE;
      break;
    case BARISTA_SCORE:
      derived_score_type = BARISTA_QVALUE;
      break;
    case BARISTA_QVALUE:
      derived_score_type = BARISTA_PEPTIDE_QVALUE;
      break;
    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case QRANKER_PEP:
    case BARISTA_PEPTIDE_QVALUE:
    case BARISTA_PEP:
    case DECOY_XCORR_PEP:
    case LOGP_WEIBULL_PEP:
    case PERCOLATOR_PEP:
    case NUMBER_SCORER_TYPES:
    case INVALID_SCORER_TYPE:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }

    match->setScore(derived_score_type, qvalue);
    scored_type_[derived_score_type] = true;

  }
  delete match_iterator;
}

/**
 * Given a hash table that maps from a score to its q-value, assign
 * q-values to all of the matches in a given collection.
 */
void MatchCollection::assignQValues(
  const map<FLOAT_T, FLOAT_T>* score_to_qvalue_hash,
  SCORER_TYPE_T score_type,
  SCORER_TYPE_T derived_score_type
){

  // Iterate over the matches filling in the q-values
  MatchIterator* match_iterator = 
    new MatchIterator(this, score_type, false);

  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    FLOAT_T score = match->getScore(score_type);

    // Retrieve the corresponding q-value.
    map<FLOAT_T, FLOAT_T>::const_iterator map_position 
      = score_to_qvalue_hash->find(score);
    if (map_position == score_to_qvalue_hash->end()) {
      carp(CARP_FATAL,
           "Cannot find q-value corresponding to score of %g.",
           score);
    }
    FLOAT_T qvalue = map_position->second;

    match->setScore(derived_score_type, qvalue);
  }
  scored_type_[derived_score_type] = true;
  delete match_iterator;
}
/**
 * Given a hash table that maps from a score to its PEP, assign
 * PEPs to all of the matches in a given collection.
 */
void MatchCollection::assignPEPs(
    const map<FLOAT_T, FLOAT_T>* score_to_pep_hash,
    SCORER_TYPE_T score_type )
{
  // Iterate over the matches filling in the q-values
  MatchIterator* match_iterator = 
    new MatchIterator(this, score_type, false);

  while(match_iterator->hasNext()){
    Match* match = match_iterator->next();
    FLOAT_T score = match->getScore(score_type);

    // Retrieve the corresponding PEP.
    map<FLOAT_T, FLOAT_T>::const_iterator map_position 
      = score_to_pep_hash->find(score);
    if (map_position == score_to_pep_hash->end()) {
      carp(CARP_FATAL,
           "Cannot find q-value corresponding to score of %g.",
           score);
    }
    FLOAT_T qvalue = map_position->second;

    /* If we're given a base score, then store the q-value.  If we're
       given a q-value, then store the peptide-level q-value. */
    SCORER_TYPE_T derived_score_type = INVALID_SCORER_TYPE;
    switch (score_type) {
    case XCORR:
      derived_score_type = DECOY_XCORR_PEP;
      break;
    case EVALUE:
      derived_score_type = DECOY_EVALUE_PEP;
      break;
    case LOGP_BONF_WEIBULL_XCORR: 
      derived_score_type = LOGP_WEIBULL_PEP;
      break;
    case PERCOLATOR_SCORE:
      derived_score_type = PERCOLATOR_PEP;
      break;
    case QRANKER_SCORE:
      derived_score_type = QRANKER_PEP;
      break;
    case BARISTA_SCORE:
      derived_score_type = BARISTA_PEP;
      break;
    // Should never reach this point.
    case SP: 
    case LOGP_WEIBULL_XCORR: 
    case DECOY_XCORR_PEPTIDE_QVALUE:
    case DECOY_XCORR_QVALUE:
    case LOGP_QVALUE_WEIBULL_XCORR:
    case LOGP_PEPTIDE_QVALUE_WEIBULL:
    case PERCOLATOR_PEPTIDE_QVALUE:
    case QRANKER_PEPTIDE_QVALUE:
    case QRANKER_PEP:
    case QRANKER_QVALUE:
    case BARISTA_PEPTIDE_QVALUE:
    case BARISTA_PEP:
    case BARISTA_QVALUE:
    case DECOY_XCORR_PEP:
    case LOGP_WEIBULL_PEP:
    case PERCOLATOR_QVALUE:
    case PERCOLATOR_PEP:
    case NUMBER_SCORER_TYPES:
    case INVALID_SCORER_TYPE:
      carp(CARP_FATAL, "Something is terribly wrong!");
    }

    match->setScore(derived_score_type, qvalue);
    scored_type_[derived_score_type] = true;

  }
  delete match_iterator;

}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


