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
#include "MatchIterator.h"
#include <string>
#include "io/MatchFileReader.h"
#include "io/SQTReader.h"
#include "util/AminoAcidUtil.h"
#include "util/Params.h"
#include "util/GlobalParams.h"
#include "util/StringUtils.h"
#include "util/WinCrux.h"

using namespace std;
using namespace Crux;

/**
 * \returns An (empty) match_collection object.
 */
void MatchCollection::init() {
  try {
    match_.reserve(10 * MILLION);
  } catch (std::bad_alloc& ba) {
    carp(CARP_DEBUG, "Bad alloc in reserve: %s", ba.what());
  }
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
  eta_ = 0;
  beta_ = 0;
  shift_ = 0;
  correlation_ = 0;
  xcorrs_.reserve(10 * MILLION);

  post_process_collection_ = false;
  top_scoring_sp_ = NULL;
  exact_pval_search_ = false;
  has_distinct_matches_ = false;
  has_decoy_indexes_ = false;
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
 * Sort the match collection by score type.
 */
void MatchCollection::sort(
  SCORER_TYPE_T score_type ///< the score type to sort by -in
  )
{
  carp(CARP_DETAILED_DEBUG, "Sorting match collection.");

  // check if we are allowed to alter match_collection
  if (iterator_lock_) {
    carp(CARP_FATAL,
         "Cannot sort a match collection when a match iterator is already instantiated");
  }

  // The order here follows that of the definition in ../objects.h.
  SCORER_TYPE_T sort_by = score_type;
  bool smaller_is_better = false;
  switch (score_type) {
  case SP:
  case XCORR:
    smaller_is_better = false;
    break;
  case EVALUE:
    smaller_is_better = true;
    break;

  case DECOY_XCORR_QVALUE:
  case DECOY_XCORR_PEPTIDE_QVALUE:
  case DECOY_XCORR_PEP:
    smaller_is_better = false;
    sort_by = XCORR;
    break;

  case DECOY_EVALUE_QVALUE:
  case DECOY_EVALUE_PEPTIDE_QVALUE:
  case DECOY_EVALUE_PEP:
    smaller_is_better = true;
    sort_by = EVALUE;
    break;

  // N.B. These are actually NEGATIVE log p-values.
  case LOGP_WEIBULL_XCORR:
  case LOGP_BONF_WEIBULL_XCORR:
  case LOGP_QVALUE_WEIBULL_XCORR:
  case LOGP_WEIBULL_PEP:
  case LOGP_PEPTIDE_QVALUE_WEIBULL:
    smaller_is_better = false;
    sort_by = LOGP_BONF_WEIBULL_XCORR;
    break;

  case PERCOLATOR_SCORE:
  case PERCOLATOR_QVALUE:
  case PERCOLATOR_PEPTIDE_QVALUE:
  case PERCOLATOR_PEP:
    smaller_is_better = false;
    sort_by = PERCOLATOR_SCORE;
    break;

  case QRANKER_SCORE:
  case QRANKER_QVALUE:
  case QRANKER_PEPTIDE_QVALUE:
  case QRANKER_PEP:
    smaller_is_better = false;
    sort_by = QRANKER_SCORE;
    break;

  case BARISTA_SCORE:
  case BARISTA_QVALUE:
  case BARISTA_PEPTIDE_QVALUE:
  case BARISTA_PEP:
    smaller_is_better = false;
    sort_by = BARISTA_SCORE;
    break;

  case DELTA_CN:
  case DELTA_LCN:
  case BY_IONS_MATCHED:
  case BY_IONS_TOTAL:
    carp(CARP_FATAL, "Cannot sort by score type %s.",
	 scorer_type_to_string(score_type));
    break;

  case TIDE_SEARCH_EXACT_PVAL:
    smaller_is_better = true;
    break;
  case TIDE_SEARCH_REFACTORED_XCORR:
    smaller_is_better = false;
    break;

  case RESIDUE_EVIDENCE_SCORE:
    smaller_is_better = false;
    break;
  case RESIDUE_EVIDENCE_PVAL:
    smaller_is_better = true;
    break;
  case BOTH_PVALUE:
    smaller_is_better = true;
    break;

  case SIDAK_ADJUSTED:
  case TIDE_SEARCH_EXACT_SMOOTHED:
    smaller_is_better = true;
    break;

  case XCORR_FIRST:
  case XCORR_SECOND:
    carp(CARP_FATAL, "Cannot sort by score type %s.",
	 scorer_type_to_string(score_type));
    
  case QVALUE_TDC:
  case QVALUE_MIXMAX:
    smaller_is_better = true;
    break;

  case NUMBER_SCORER_TYPES:
  case INVALID_SCORER_TYPE:
    carp(CARP_FATAL, "Cannot sort by score type %s.",
	 scorer_type_to_string(score_type));
  }

  // Don't sort if it's already sorted.
  if (last_sorted_ == sort_by) {
    return;
  }

  if (!scored_type_[sort_by]) {
    carp(CARP_WARNING, "Cannot sort MatchCollection (does not have %s scores)",
	 scorer_type_to_string(sort_by));
    return;
  }

  if (smaller_is_better) {
    carp(CARP_DEBUG, "Sorting in ascending order by %s.",
	 scorer_type_to_string(score_type));
  } else {
    carp(CARP_DEBUG, "Sorting in descending order by %s.",
	 scorer_type_to_string(score_type));
  }	 

  // Do the sort.
  Match::ScoreComparer comparer(sort_by, smaller_is_better);
  std::sort(match_.begin(), match_.end(), comparer);
  last_sorted_ = sort_by;
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
  if (last_sorted_ != score_type) {
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
  const string& file_path,  ///< File path to set
  bool overwrite ///< Overwrite existing files?
  ) {
  int file_idx = -1;
  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    if (overwrite || (*i)->getFileIndex() == -1) {
      file_idx = (*i)->setFilePath(file_path);
    }
  }
  return file_idx;
}

/**
 *\returns true, if there is a  match_iterators instantiated by match collection 
 */
bool MatchCollection::getIteratorLock() {
  return iterator_lock_;
}

/**
 *\returns the total match objects avaliable in current match_collection
 */
int MatchCollection::getMatchTotal() {
  return match_.size();
}

bool MatchCollection::getHasDistinctMatches() {
  return has_distinct_matches_;
}

void MatchCollection::setHasDistinctMatches(bool distinct) {
  has_distinct_matches_ = distinct;
}

bool MatchCollection::hasDecoyIndexes() const {
  return has_decoy_indexes_;
}

void MatchCollection::setHasDecoyIndexes(bool value) {
  has_decoy_indexes_ = value;
}

void MatchCollection::setExperimentSize(int size) {
  experiment_size_ = size;
}

/**
 * \returns The total number of peptides searched for this spectrum,
 * target peptides for a target collection or decoy peptides for a
 * decoy collection.
 */
int MatchCollection::getExperimentSize() {
  return experiment_size_;
}

/**
 * Sets the total number of target peptides searched for this
 * spectrum.  Only needs to be used by decoy match collections.
 */
void MatchCollection::setTargetExperimentSize(int numMatches) {
  target_experiment_size_ = numMatches;
}

/**
 * \returns The number of target matches that this spectrum had.
 * Different than getExperimentSize() for decoy match collections.
 */
int MatchCollection::getTargetExperimentSize() {
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
  FILE* output,
  const string& ms2file
  ){
  if (output == NULL) {
    return;
  }
  time_t hold_time;
  ENZYME_T enzyme = GlobalParams::getEnzyme();
  const char* enz_str = enzyme_type_to_string(enzyme);
  string database = Params::GetString("protein-database");

  MASS_TYPE_T isotopic_mass_type = GlobalParams::getIsotopicMass();
  MASS_TYPE_T fragment_mass_type = GlobalParams::getFragmentMass();

  const char* isotopic_mass;
  const char* fragment_mass;
  DIGEST_T digest = GlobalParams::getDigestion();
  int max_num_internal_cleavages;
  int min_number_termini;
  int missed_cleavage = GlobalParams::getMissedCleavages();
  if (missed_cleavage){
    max_num_internal_cleavages = GlobalParams::getMaxLength();
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
          StringUtils::RTrim(string(ctime(&hold_time))).c_str(),
          xmlns, xsi, xmlns, schema_location);
  fprintf(output, "<msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" "
          "msModel=\"%s\" msIonization=\"%s\" msAnalyzer=\"%s\" "
          "msDectector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\" >\n",
          "NA", // TODO, dummy value
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
          "NA", // TODO, dummy value
          "Crux",
          isotopic_mass, // isotopic mass type is precursor mass type?
          fragment_mass,
          "NA", // TODO, dummy value
          "NA",
          1 // TODO, dummy value
          );
  fprintf(output, "<search_database local_path=\"%s\" type=\"%s\" />\n", 
          absolute_database_path, 
          "AA");
  fprintf(output, "<enzymatic_search_constraint enzyme=\"%s\" "
          "max_num_internal_cleavages=\"%i\" min_number_termini=\"%i\"/>\n",
          enz_str,
          max_num_internal_cleavages,
          min_number_termini);
#ifndef DARWIN
  free(absolute_database_path);
#endif


  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A'+ ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = GlobalParams::getIsotopicMass();
  int aa = 0;

  // static amino acid modifications
  for (char aa = 'A'; aa <= 'Z'; aa++) {
    vector<const ModificationDefinition*> staticMods = ModificationDefinition::StaticMods(aa);
    double aaMass = AminoAcidUtil::GetMass(aa, isotopic_type == MONO);
    for (vector<const ModificationDefinition*>::const_iterator i = staticMods.begin();
         i != staticMods.end();
         i++) {
      double modMass = (*i)->DeltaMass();
      double totalMass = aaMass + modMass;
      string termString;
      switch ((*i)->Position()) {
        case PEPTIDE_N: termString = "peptide_terminus=\"n\" "; break;
        case PEPTIDE_C: termString = "peptide_terminus=\"c\" "; break;
      }
      fprintf(output, "<aminoacid_modification aminoacid=\"%c\" mass=\"%s\" "
                      "massdiff=\"%s\" variable=\"N\" %s/>\n",
              aa,
              StringUtils::ToString(totalMass, Params::GetInt("mass-precision")).c_str(),
              StringUtils::ToString(modMass, Params::GetInt("mod-precision")).c_str(),
              termString.c_str());
    }
  }

  // variable amino acid modifications
  vector<const ModificationDefinition*> varMods = ModificationDefinition::VarMods();
  for (vector<const ModificationDefinition*>::const_iterator i = varMods.begin();
       i != varMods.end();
       i++) {
    for (set<char>::const_iterator j = (*i)->AminoAcids().begin();
         j != (*i)->AminoAcids().end();
         j++) {
      double modMass = (*i)->DeltaMass();
      double aaMass = AminoAcidUtil::GetMass(*j, isotopic_type == MONO);
      double totalMass = aaMass + modMass;
      string termString;
      switch ((*i)->Position()) {
        case PEPTIDE_N: termString = "peptide_terminus=\"n\" "; break;
        case PEPTIDE_C: termString = "peptide_terminus=\"c\" "; break;
      }
      fprintf(output, "<aminoacid_modification aminoacid=\"%c\" mass=\"%s\" "
                      "massdiff=\"%s\" variable=\"Y\" %s/>\n",
              *j,
              StringUtils::ToString(totalMass, Params::GetInt("mass-precision")).c_str(),
              StringUtils::ToString(modMass, Params::GetInt("mod-precision")).c_str(),
              termString.c_str());
    }
  }

  for (vector<const Param*>::const_iterator i = Params::Begin(); i != Params::End(); i++) {
    string name = (*i)->GetName();
    if ((*i)->IsVisible() &&
        name != "mod" && name != "cmod" && name != "nmod") {
      fprintf(output, "<parameter name=\"%s\" value=\"%s\"/>\n",
              name.c_str(), (*i)->GetString().c_str());
    }
  }
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

  MASS_TYPE_T mass_type = GlobalParams::getIsotopicMass();
  char temp_str[64];
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tPrecursorMasses\t%s\n", temp_str);
  
  mass_type = get_mass_type_parameter("fragment-mass");
  mass_type_to_string(mass_type, temp_str);
  fprintf(output, "H\tFragmentMasses\t%s\n", temp_str); //?????????

  double tol = Params::GetDouble("precursor-window");
  fprintf(output, "H\tAlg-PreMasTol\t%.1f\n",tol);
  fprintf(output, "H\tAlg-FragMassTol\t%.2f\n", 
          Params::GetDouble("mz-bin-width") / 2.0);
  fprintf(output, "H\tAlg-XCorrMode\t0\n");

  fprintf(output, "H\tComment\tpreliminary algorithm %s\n", 
          Params::GetString("prelim-score-type").c_str());

  fprintf(output, "H\tComment\tfinal algorithm %s\n",
          Params::GetString("score-type").c_str());

  int aa = 0;
  char aa_str[2];
  aa_str[1] = '\0';
  int alphabet_size = (int)'A' + ((int)'Z'-(int)'A');
  MASS_TYPE_T isotopic_type = GlobalParams::getIsotopicMass();

  for(aa = (int)'A'; aa < alphabet_size -1; aa++){
    aa_str[0] = (char)aa;
    double mod = Params::GetDouble(aa_str);
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
    string aa_list_str = aamod->getAAListString();
    char aa_symbol = aamod->getSymbol();
    double mass_dif = aamod->getMassChange();

    fprintf(output, "H\tDiffMod\t%s%c=%+.2f\n", aa_list_str.c_str(), 
            aa_symbol, mass_dif);
  }
  num_mods = get_c_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aamod->getSymbol();

    fprintf(output, "H\tComment\tMod %c is a C-terminal modification\n",
            aa_symbol);
  }

  num_mods = get_n_mod_list(&aa_mod_list);
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    AA_MOD_T* aamod = aa_mod_list[mod_idx];
    char aa_symbol = aamod->getSymbol();

    fprintf(output, "H\tComment\tMod %c is a N-terminal modification\n",
            aa_symbol);
  }



  //for letters in alphabet
  //  double mod = Params::GetDouble(letter);
  //  if mod != 0
  //     double mass = mod + getmass(letter);
  //     fprintf(output, "H\tStaticMod\t%s=%.3f\n", letter, mass);
  //  fprintf(output, "H\tStaticMod\tC=160.139\n");
  fprintf(output, "H\tAlg-DisplayTop\t%d\n", 
          //          Params::GetInt("max-sqt-result")); 
          Params::GetInt("top-match")); 
  // this is not correct for an sqt from analzyed matches

  ENZYME_T enzyme = get_enzyme_type_parameter("enzyme");
  DIGEST_T digestion = get_digest_type_parameter("digestion");
  const char* enz_str = enzyme_type_to_string(enzyme);
  const char* dig_str = digest_type_to_string(digestion);
  char custom_str[SMALL_BUFFER];
  if( enzyme == CUSTOM_ENZYME){
    string rule = Params::GetString("custom-enzyme");
    sprintf(custom_str, ", custom pattern: %s", rule.c_str());
  }else{
    custom_str[0] = 0;
  }
  fprintf(output, "H\tEnzymeSpec\t%s-%s%s\n", enz_str, dig_str, custom_str);
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
void MatchCollection::printMultiSpectraXml(PepXMLWriter* output) {
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
    if (!is_decoy) {
      int* ranks = new int[NUMBER_SCORER_TYPES];
      ranks[XCORR] = -1;
      char* peptide_sequence = cur_match->getSequence();
      char* mod_peptide_sequence = cur_match->getModSequenceStrWithMasses(
                             get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = cur_match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, protein_descriptions);
      for (int score_idx = 0; score_idx < NUMBER_SCORER_TYPES; score_idx++) {
        if (scored_type_[score_idx]) {
          scores[score_idx] = cur_match->getScore((SCORER_TYPE_T)score_idx);
          ranks[score_idx] = cur_match->getRank((SCORER_TYPE_T)score_idx);
        }
      }
      unsigned num_matches = (!isDecoy()) ? getTargetExperimentSize() : getExperimentSize();
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        cur_match->getNeutralMass(),
        cur_match->getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->calcModifiedMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        scored_type_,
        scores,
        num_matches);
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
  if (scored_type_[SP]) {
    scores_computed[SP] = true;
  }
  if (scored_type_[DELTA_CN]) {
    scores_computed[DELTA_CN] = true;
  }
  if (scored_type_[BY_IONS_MATCHED]) {
    scores_computed[BY_IONS_MATCHED] = true;
  }
  if (scored_type_[BY_IONS_TOTAL]) {
    scores_computed[BY_IONS_TOTAL] = true;
  }
  scores_computed[TIDE_SEARCH_EXACT_PVAL] = exact_pval_search_;
  scores_computed[TIDE_SEARCH_REFACTORED_XCORR] = exact_pval_search_;
  scores_computed[main_score] = !exact_pval_search_;

  double* scores = new double[NUMBER_SCORER_TYPES];
  int* ranks = new int[NUMBER_SCORER_TYPES];

  Match* match = NULL;
  // create match iterator
  // true: return match in sorted order of main_score type
  MatchIterator* match_iterator = new MatchIterator(this, main_score, true);
  // iterate over matches
  while (match_iterator->hasNext()) {
    match = match_iterator->next();
    int cur_rank = match->getRank(main_score);   
    if (scored_type_[XCORR]) {
      ranks[XCORR] = match->getRank(XCORR);
    }
    if (scored_type_[SP]) {
      ranks[SP] = match->getRank(SP);
      scores[SP] = match->getScore(SP);
    }
    if (scored_type_[DELTA_CN]) {
      scores[DELTA_CN] = match->getScore(DELTA_CN);
    }
    if (scored_type_[BY_IONS_MATCHED]) {
      scores[BY_IONS_MATCHED] = match->getScore(BY_IONS_MATCHED);
    }
    if (scored_type_[BY_IONS_TOTAL]) {
      scores[BY_IONS_TOTAL] = match->getScore(BY_IONS_TOTAL);
    }
    // print if we haven't reached the limit
    // or if we are at the limit but this match is a tie with the last
    if (count < top_match || last_rank == cur_rank) {
      char* peptide_sequence = match->getSequence();
      char* mod_peptide_sequence = match->getModSequenceStrWithMasses(
                            get_mass_format_type_parameter("mod-mass-format"));
      Peptide* peptide = match->getPeptide();
      char* flanking_aas = peptide->getFlankingAAs();
      int num_proteins = peptide->getProteinInfo(protein_ids, 
                                                 protein_descriptions);
      for(int score_idx=0; score_idx < NUMBER_SCORER_TYPES; score_idx++){
        if(scored_type_[score_idx]) {
          scores[score_idx] = match->getScore((SCORER_TYPE_T)score_idx);
        }
      }   
      unsigned num_matches = (!isDecoy()) ? getTargetExperimentSize() : getExperimentSize();
      output->writePSM(spectrum->getFirstScan(),
        spectrum->getFilename(),
        zstate_.getNeutralMass(),
        zstate_.getCharge(),
        ranks,
        peptide_sequence,
        mod_peptide_sequence,
        peptide->calcModifiedMass(),
        num_proteins,
        flanking_aas,
        protein_ids,
        protein_descriptions,
        scores_computed,
        scores,
        num_matches);
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
  delete [] scores_computed;
  delete [] scores;
  delete [] ranks;

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

      match->printTab(this, output, spectrum, num_target_matches, num_decoy_matches);
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
) {
  carp(CARP_DETAILED_DEBUG, "Writing matches to file");

  // if file location is target (i.e. tdc=T), print all to target
  MatchFileWriter* decoy_file = decoy_tab_file;
  if (Params::GetString("estimation-method") == "tdc") {
    decoy_file = tab_file;
  }

  // for each match, get spectrum info, determine if decoy, print
  for (vector<Crux::Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    Match* cur_match = *i;
    bool is_decoy = cur_match->getNullPeptide();
    FLOAT_T num_psm_per_spec = cur_match->getLnExperimentSize();
    num_psm_per_spec = expf(num_psm_per_spec) + 0.5; // round to nearest int
    int num_target_psm_per_spec = cur_match->getTargetExperimentSize();

    cur_match->printTab(this, (is_decoy ? decoy_file : tab_file),
                        cur_match->getSpectrum(),
                        (int)num_psm_per_spec, num_target_psm_per_spec);
  }
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
    carp(CARP_DEBUG, "parsed PSM: %d", match_.size());
  }
  
  return true;
}

/**
 * \brief Calculate the delta_cn of each match and populate the field.
 * 
 * Delta_cn is the normalized difference between xcorrs of different ranks.
 * match[i] = (match[i] - match[i+1]) / match[i].
 * Sorts match_collection by xcorr, if necessary.
 * 
 */
bool MatchCollection::calculateDeltaCn() {
  if (scored_type_[DELTA_CN] && scored_type_[DELTA_LCN]) {
    return false;
  }

  SCORER_TYPE_T type = XCORR;
  if (!scored_type_[type]) {
    type = TIDE_SEARCH_EXACT_PVAL;
    if (!scored_type_[type]) {
      carp(CARP_WARNING,
           "Delta_cn not calculated because match collection not scored for xcorr");
      return false;
    }
  }

  vector<FLOAT_T> scores;
  for (vector<Match*>::iterator i = match_.begin(); i != match_.end(); i++) {
    scores.push_back((*i)->getScore(type));
  }
  vector< pair<FLOAT_T, FLOAT_T> > deltaCns = calculateDeltaCns(scores, type);

  for (size_t i = 0; i < deltaCns.size(); i++) {
    match_[i]->setScore(DELTA_CN, deltaCns[i].first);
    match_[i]->setScore(DELTA_LCN, deltaCns[i].second);
  }

  scored_type_[DELTA_CN] = scored_type_[DELTA_LCN] = true;
  return true;
}

vector< pair<FLOAT_T, FLOAT_T> > MatchCollection::calculateDeltaCns(
  vector<FLOAT_T> scores,
  SCORER_TYPE_T type
) {
  vector< pair<FLOAT_T, FLOAT_T> > deltaCns(scores.size(), make_pair(0, 0));
  if (scores.empty()) {
    return deltaCns;
  }

  map< const FLOAT_T*, pair<FLOAT_T, FLOAT_T>* > scorePtrs;
  for (size_t i = 0; i < scores.size(); i++) {
    scorePtrs[&scores[i]] = &deltaCns[i];
  }

  if (type != TIDE_SEARCH_EXACT_PVAL) {
    // Higher is better - sort descending
    std::sort(scores.begin(), scores.end(), std::greater<FLOAT_T>());
  } else {
    // Lower is better - sort ascending
    std::sort(scores.begin(), scores.end(), std::less<FLOAT_T>());
  }
  
  int top_match = Params::GetInt("top-match");

//  FLOAT_T last = scores.back();
  FLOAT_T last = scores[min(top_match, (int) scores.size()) - 1];
  for (vector<FLOAT_T>::const_iterator i = scores.begin(); i != scores.end(); i++) {
    vector<FLOAT_T>::const_iterator next = (i != scores.end() - 1) ? i + 1 : i;
    FLOAT_T deltaCn, deltaLCn;
    switch (type) {
      default:
        deltaCn = (*i - *next) / max(*i, (FLOAT_T)1.0);
        deltaLCn = (*i - last) / max(*i, (FLOAT_T)1.0);
        break;
      case TIDE_SEARCH_EXACT_PVAL:
        deltaCn = -log10(*i) + log10(*next);
        deltaLCn = -log10(*i) + log10(last);
        break;
    }
    scorePtrs[&*i]->first = deltaCn;
    scorePtrs[&*i]->second = deltaLCn;
  }
  return deltaCns;
}

/**********************************
 * match_collection get, set methods
 **********************************/

/**
 * \returns true if the match_collection only contains decoy matches,
 * else (all target or mixed) returns false.
 */
bool MatchCollection::isDecoy() {
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

/**
 * Extract a given type of score into a vector.
 */
vector<FLOAT_T> MatchCollection::extractScores(
  SCORER_TYPE_T score_type ///< Type of score to extract.
) const {
  size_t numMatches = match_.size();
  vector<FLOAT_T> scores(numMatches, 0);
  for (size_t i = 0; i < numMatches; i++) {
    scores[i] = match_[i]->getScore(score_type);
  }
  return scores;
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

    FLOAT_T qvalue;
    // If the score is not a number, punt.
    if ( isinf(score) || isnan(score) ) {
      carp(CARP_DEBUG, "Found inf or nan score.");
      qvalue = numeric_limits<double>::quiet_NaN();
    } else {
      // Retrieve the corresponding q-value.
      map<FLOAT_T, FLOAT_T>::const_iterator map_position 
	= score_to_qvalue_hash->find(score);
      if (map_position == score_to_qvalue_hash->end()) {
	carp(CARP_FATAL,
	     "Cannot find q-value corresponding to score of %g.",
	     score);
      }
      qvalue = map_position->second;
    }
    match->setScore(derived_score_type, qvalue);
  }
  scored_type_[derived_score_type] = true;
  delete match_iterator;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


