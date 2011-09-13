#include "SpectralCounts.h"
#include "crux-utils.h"
#include "OutputFiles.h"
#include "ProteinPeptideIterator.h"
#include "SpectrumCollectionFactory.h"
#include "MatchCollectionIterator.h"

using namespace std;

/**
 * Default constructor.
 */
SpectralCounts::SpectralCounts() 
  : output_(NULL), 
    threshold_(0),
    unique_mapping_(false),
    quantitation_(PROTEIN_QUANT_LEVEL),
    parsimony_(PARSIMONY_NONE),
    measure_(MEASURE_SIN),
    bin_width_(0),
    peptide_scores_(Peptide::lessThan),
    protein_scores_(protein_id_less_than),
    protein_supporting_peptides_(protein_id_less_than),
    protein_meta_protein_(protein_id_less_than),
    meta_mapping_(comparePeptideSets),
    meta_protein_scores_(compareMetaProteins),
    meta_protein_ranks_(compareMetaProteins){
}

/**
 * Default destructor.
 */
SpectralCounts::~SpectralCounts() {
  delete output_;
}



/**
 * Given a collection of scored PSMs, print a list of proteins
 * ranked by their a specified score. Spectral-counts supports two
 * types of quantification: Normalized Spectral Abundance Factor (NSAF)
 * and Normalized Spectral Index (SIN). 
 * \returns 0 on successful completion.
 */
int SpectralCounts::main(int argc, char** argv) {

  const char* option_list[] = {
    "verbosity",
    "parameter-file",
    "parsimony",
    "threshold",
    "input-ms2",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "quant-level",
    "measure"
  };
  const char* argument_list[] = {
    "input PSMs",
    "protein database"
  };

  int num_options = sizeof(option_list) / sizeof(char*);
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize(argument_list, num_arguments,
             option_list, num_options, argc, argv);

  getParameterValues(); // all the get_<type>_parameter calls here

  // open output files
  output_ = new OutputFiles(this);
  output_->writeHeaders();

  // get a set of matches that pass the threshold
  filterMatches();
  carp(CARP_INFO, "Number of matches passed the threshold %i", 
       matches_.size());

  // get a set of peptides
  getPeptideScores();
  if( unique_mapping_ ){
    makeUniqueMapping();
  }
  carp(CARP_INFO, "Number of Unique Peptides %i", peptide_scores_.size());

  // quantify at either the peptide or protein level
  if( quantitation_ == PEPTIDE_QUANT_LEVEL ){ // peptide level
    normalizePeptideScores();
    output_->writeRankedPeptides(peptide_scores_);

  } else if( quantitation_ == PROTEIN_QUANT_LEVEL ){ // protein level
    
    getProteinScores();
    normalizeProteinScores();
    carp(CARP_INFO, "Number of Proteins %i", protein_scores_.size());
        
    if( parsimony_ != PARSIMONY_NONE ){ //if parsimony is not none
      getProteinToPeptides();
      getMetaMapping();
      getProteinToMetaProtein();
      carp(CARP_INFO, "Number of meta proteins %i", meta_mapping_.size());
      
      if( parsimony_ == PARSIMONY_GREEDY ){ //if parsimony is greedy
        performParsimonyAnalysis();
      }
      getMetaScores();
      getMetaRanks();
    }
    
    output_->writeRankedProteins(protein_scores_, meta_protein_ranks_, 
                                 protein_meta_protein_);

  } else {
    carp(CARP_FATAL, "Invalid quantification level.");
  }
  
  return 0;
}

/**
 * Collect all the parameter/option values we will need and store them
 * as member variables.
 */
void SpectralCounts::getParameterValues(){
  psm_file_ = get_string_parameter_pointer("input PSMs");
  threshold_ = get_double_parameter("threshold");
  database_name_ = get_string_parameter_pointer("protein database");
  unique_mapping_ = get_boolean_parameter("unique-mapping");
  quantitation_ = get_quant_level_type_parameter("quant-level");
  parsimony_ = get_parsimony_type_parameter("parsimony");
  measure_ = get_measure_type_parameter("measure");
  bin_width_ = get_double_parameter("mz-bin-width");
}

/**
 * For every protein that can be mapped from the set of 
 * peptides in PeptideToScore map, enter the protein and
 * the set of identified peptides it maps to, into 
 * ProteinToPeptide
 */
void SpectralCounts::getProteinToPeptides(){
  for (PeptideToScore::iterator pep_it = peptide_scores_.begin();
       pep_it != peptide_scores_.end(); ++pep_it){
    Peptide* peptide = pep_it->first;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PeptideSrc* peptide_src = 
        peptide_src_iterator_next(peptide_src_iterator);
      Protein* protein = peptide_src->getParentProtein();
      if (protein_supporting_peptides_.find(protein) == 
          protein_supporting_peptides_.end()){
        PeptideSet newset(Peptide::lessThan);
        protein_supporting_peptides_.insert(make_pair(protein, newset));
      }
      protein_supporting_peptides_[protein].insert(peptide);
    }
    free(peptide_src_iterator);
  }
}


/**
 * Enters the mapping of protein to its metaProtein
 * into ProteinToMetaProtein. MetaProteins are retreieved
 * from MetaMapping
 *
 */
void SpectralCounts::getProteinToMetaProtein(){
  // for every meta protein
  for (MetaMapping::iterator meta_protein_it = meta_mapping_.begin();
       meta_protein_it != meta_mapping_.end(); ++meta_protein_it){
    MetaProtein proteins = meta_protein_it->second;
    // for every protein in the meta protein
    for (MetaProtein::iterator proteins_it = proteins.begin();
         proteins_it != proteins.end(); ++proteins_it){
      // create a mapping of protein to meta protein
      protein_meta_protein_.insert(make_pair((*proteins_it), proteins));
    }
  }
}

/**
 * A score for each protein is calculated by summing 
 * the scores of each peptide that belongs to a protein
 */
void SpectralCounts::getProteinScores(){

  // iterate through each peptide
  for (PeptideToScore::iterator pep_it = peptide_scores_.begin();
       pep_it != peptide_scores_.end(); ++pep_it){
    Peptide* peptide = pep_it->first;
    FLOAT_T pep_score = pep_it->second;
    PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator =
      new_peptide_src_iterator(peptide);
    while( peptide_src_iterator_has_next(peptide_src_iterator)) {
      PeptideSrc* peptide_src = 
        peptide_src_iterator_next(peptide_src_iterator);
      Protein* protein = peptide_src->getParentProtein();
      if (protein_scores_.find(protein) == protein_scores_.end()){
        protein_scores_.insert(make_pair(protein, 0.0));
      }
      protein_scores_[protein] += pep_score;
    }
    free(peptide_src_iterator);
  }
}

/**
 * Takes the PeptideToScores map and updates all
 * values with normalized values. Normalized by sum
 * of all scores and then by the peptide length
 *
 */
void SpectralCounts::normalizePeptideScores()
{
  carp(CARP_INFO, "Normalizing peptide scores");
  FLOAT_T total = 0.0;

  // calculate sum of all scores
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it){
    FLOAT_T score = it->second;
    total += score;
  }

  // normalize by sum of scores and length
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it){
    FLOAT_T score = it->second;
    Peptide* peptide = it->first;
    it->second = score / total / peptide->getLength();

  }

}

/**
 * Changes the scores in protien_scores_ to either be divided by the
 * sum of all scores times the peptide length (SIN, NSAF) or to be the
 * final emPAI score.
 */
void SpectralCounts::normalizeProteinScores(){
  if( measure_ == MEASURE_EMPAI ){
    computeEmpai();
  } else {

    carp(CARP_INFO, "Normalizing protein scores");
    FLOAT_T total = 0.0;
    
    // calculate sum of all scores
    for (ProteinToScore::iterator it = protein_scores_.begin();
         it != protein_scores_.end(); ++it){
      FLOAT_T score = it->second;
      total += score;
    }
    
    // normalize by sum of all scores and length
    for (ProteinToScore::iterator it = protein_scores_.begin();
         it != protein_scores_.end(); ++it){
      FLOAT_T score = it->second;
      Protein* protein = it->first;
      it->second = score / total / protein->getLength();
    }
  }
}

/**
 * Computes the 10^(observed/total) - 1 score for each protein.  Assumes
 * that protein_scores_ has been populated with the count of observed
 * unique peptides for each protein.
 */
void SpectralCounts::computeEmpai(){
  PeptideConstraint* constraint = PeptideConstraint::newFromParameters();

  ProteinToScore::iterator it = protein_scores_.begin();;
  for (; it != protein_scores_.end(); ++it){
    FLOAT_T observed_peptides = it->second;
    Protein* protein = it->first;
    ProteinPeptideIterator* iter = new ProteinPeptideIterator(protein, 
                                                              constraint);
    FLOAT_T possible_peptides = iter->getTotalPeptides();

    it->second = pow(10, (observed_peptides / possible_peptides)) - 1;

    delete iter;
  }
  PeptideConstraint::free(constraint);
}

/**
 * For the spectrum associated with the match, sum the intensities of
 * all b and y ions that are not modified.
 * \return The sum of unmodified b and y ions.
 */
FLOAT_T SpectralCounts::sumMatchIntensity(Match* match,
                                        SpectrumCollection* spectra)
{
  FLOAT_T match_intensity = 0;
  char* peptide_seq = match->getSequence();
  MODIFIED_AA_T* modified_sequence = match->getModSequence();
  int charge = match->getCharge();
  Spectrum* temp = match->getSpectrum();
  int scan = temp->getFirstScan();
  Spectrum* spectrum = spectra->getSpectrum(scan);
  Ion* ion;
  SCORER_TYPE_T score_type = XCORR;
  IonConstraint* ion_constraint =
    IonConstraint::newIonConstraintSmart(score_type, charge);
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);
  ion_series->update(peptide_seq, modified_sequence);
  ion_series->predictIons();
  for (IonIterator ion_it = ion_series->begin();
       ion_it != ion_series->end(); ++ion_it){
    ion = (*ion_it);
    if (ion -> getType() == B_ION || ion -> getType() == Y_ION){
      if (!ion->isModified()){
        Peak * peak = spectrum->getNearestPeak(ion->getMassZ(),
                                                bin_width_);
        if (peak != NULL){
          match_intensity += peak->getIntensity();
        }
      }
    }
  }
  delete ion_series;
  free(peptide_seq);

  return match_intensity;
}


/**
 * Generate a score for each peptide in the set of matches.  Populate
 * peptide_scores_ which becomes a unique set of peptides, each with a
 * score.
 *
 * For SIN the score is the sum of intensites of b and y ions (without
 * H2O modifications).  Intensites are taken from the .ms2 file.
 *
 * For NSAF, the score is the number of matches (spectra) for each peptide.
 * For EMPAI, the score is 1 as it only reqires the number of peptides
 * observed per protein.
 */
void SpectralCounts::getPeptideScores()
{
  SpectrumCollection* spectra = NULL;

  // for SIN, parse out spectrum collection from ms2 fiel
  if( measure_ == MEASURE_SIN ){
    spectra = SpectrumCollectionFactory::create(get_string_parameter_pointer("input-ms2"));
  }

  for(set<Match*>::iterator match_it = matches_.begin();
      match_it != matches_.end(); ++match_it){

    FLOAT_T match_intensity = 1; // for NSAF just count each for the peptide/

    Match* match = (*match_it);
    // for sin, calculate total ion intensity for match by
    // summing up peak intensities
    if( measure_ == MEASURE_SIN ){
      match_intensity = sumMatchIntensity(match, spectra);
    }

    // add ion_intensity to peptide scores
    Peptide* peptide = match->getPeptide();
    if (peptide_scores_.find(peptide) ==  peptide_scores_.end()){
      peptide_scores_.insert(make_pair(peptide, 0.0));
    }
    peptide_scores_[peptide] += match_intensity;
  }

  if( measure_ == MEASURE_SIN ){
    delete spectra;
  }

  // for emPAI we just need a count of unique peptides
  if( measure_ == MEASURE_EMPAI ){
    PeptideToScore::iterator itr = peptide_scores_.begin();
    for(; itr != peptide_scores_.end(); ++itr){
      itr->second = 1; // count each unique peptide once
    }
  }
}


/**
 * Create a set of matches, all with an XCORR rank == 1 and all of which
 * have a qvalue score lower than user-specified threshold.
 */
void SpectralCounts::filterMatches()
{
  // get input file directory
  char** path_info = parse_filename_path(psm_file_.c_str());
  if( path_info[1] == NULL ){
    path_info[1] = my_copy_string(".");
  }

  // create match collection
  int decoy_count = 0;
  MatchCollectionIterator* match_collection_it 
    = new MatchCollectionIterator(path_info[1], 
                                  database_name_.c_str(), &decoy_count);

  matches_.clear();
  MatchIterator* match_iterator = NULL;
  MatchCollection* match_collection = NULL;
  bool qualify = false;
  while (match_collection_it->hasNext()){
    
    match_collection = match_collection_it->next();
    match_iterator = new MatchIterator(match_collection, XCORR, TRUE);
    
    // figure out which qvalue we are using
    SCORER_TYPE_T qval_type = get_qval_type(match_collection);
    if( qval_type == INVALID_SCORER_TYPE ){
      carp(CARP_FATAL, "The matches in %s do not have q-values from percolator,"
           " q-ranker, or compute-q-values.\n", psm_file_.c_str());
    }

    while(match_iterator->hasNext()){
      Match* match = match_iterator->next();
      qualify = false;
      if (match->getRank(XCORR) != 1){
        continue;
      }
      // find a qvalue score lower than threshold
      if (match->getScore(qval_type) != FLT_MIN &&
          match->getScore(qval_type) <= threshold_)  {
        matches_.insert(match);
      }
    } // next match
  } // next file
  free(path_info[1]);
  free(path_info[0]);
  free(path_info);
  // cannot delete the mach_collection_iterator b/c it deletes the
  // database which deletes the proteins which are still in use after
  // filterMatches returns.
  // TODO create a database in calling function, pass to
  // filterMatches, have a MatchCollectionIterator constructor that
  // takes a database
}

/**
 * Figures out which kind of q-value was scored for this match collection.
 * \returns PERCOLATOR_QVALUE, QRANKER_QVALUE, or DECOY_XCORR_QVALUE
 * if any of those were scored or INVALID_SCORER_TYPE if none were scored. 
 */
SCORER_TYPE_T SpectralCounts::get_qval_type(
  MatchCollection* match_collection){
  SCORER_TYPE_T scored_type = INVALID_SCORER_TYPE;

  if(match_collection->getScoredType(PERCOLATOR_QVALUE)){
    scored_type =  PERCOLATOR_QVALUE;
  } else if(match_collection->getScoredType(QRANKER_QVALUE)){
    scored_type = QRANKER_QVALUE;
  } else if(match_collection->getScoredType(DECOY_XCORR_QVALUE)){
    scored_type = DECOY_XCORR_QVALUE;
  }

  return scored_type;
}

/**
 * Fills in the MetaMapping with entries of set of 
 * peptides that can be found in every protein in
 * the meta protein
 *
 */
void SpectralCounts::getMetaMapping(){
  carp(CARP_INFO, "Creating a mapping of meta protein to peptides");
  int count = 0;
  for(ProteinToPeptides::iterator prot_it= protein_supporting_peptides_.begin();
       prot_it != protein_supporting_peptides_.end(); ++prot_it){
    Protein* protein = prot_it->first;
    PeptideSet pep_set = prot_it->second;

    if (meta_mapping_.find(pep_set) == meta_mapping_.end()){
      MetaProtein meta_protein(protein_id_less_than);
      count++;
      meta_mapping_.insert(make_pair(pep_set, meta_protein));
    }
    meta_mapping_[pep_set].insert(protein);
  }

}

/**
 * Takes a mapping of set of peptides to meta proteins and 
 * a mapping of protein to scores, and finds the largest score
 * of the meta proteins for each protein. The object returned
 * is a mapping of MetaProteins to the highest score
 *
 */
void SpectralCounts::getMetaScores(){
  carp(CARP_INFO, "Finding scores of meta proteins");
  for (MetaMapping::iterator meta_it = meta_mapping_.begin();
       meta_it != meta_mapping_.end(); ++meta_it ){
    MetaProtein proteins = (*meta_it).second;
    FLOAT_T top_score = -1.0;
    for (MetaProtein::iterator protein_it = proteins.begin();
         protein_it != proteins.end(); ++protein_it){
      Protein* protein = (*protein_it);
      FLOAT_T score = protein_scores_[protein];
      top_score = max(score, top_score);
    }
    meta_protein_scores_.insert(make_pair(proteins, top_score));
  }

}

/**
 * Takes a mapping of MetaProteins to scores and returns 
 * a mapping of set of peptides to rank
 *
 */
void SpectralCounts::getMetaRanks(){
  carp(CARP_INFO, "Finding ranks of meta proteins");
  vector< pair<FLOAT_T, MetaProtein> > metaVector;
  for (MetaToScore::iterator meta_it = meta_protein_scores_.begin();
       meta_it != meta_protein_scores_.end(); ++meta_it){
    MetaProtein proteins = (*meta_it).first;
    FLOAT_T score = (*meta_it).second;
    metaVector.push_back(make_pair(score, proteins));
  }
  sort(metaVector.begin(), metaVector.end());
  reverse(metaVector.begin(), metaVector.end());

  int cur_rank = 1;
  for (vector< pair<FLOAT_T, MetaProtein> >::iterator
         vector_it = metaVector.begin();
       vector_it != metaVector.end(); ++vector_it){
    MetaProtein proteins = (*vector_it).second;
    meta_protein_ranks_.insert(make_pair(proteins, cur_rank));
    cur_rank++;
  }

}

/**
 * Greedily finds a peptide-to-protein mapping where each
 * peptide is only mapped to a single meta-protein. 
 *
 * Would of been better to implement with priority queue w/
 * adjancency lists: O(n*log(n)) but input size should be
 * small enough that performance should not be an issue.
 */
void SpectralCounts::performParsimonyAnalysis(){
  carp(CARP_INFO, "Performing Greedy Parsimony analysis");
  MetaMapping result(comparePeptideSets);
  vector< pair<PeptideSet, MetaProtein > > peps_vector;

  // get all meta mappings into a vector 
  for (MetaMapping::iterator meta_iter = meta_mapping_.begin();
       meta_iter != meta_mapping_.end(); ++meta_iter){
    peps_vector.push_back((*meta_iter));
  }


  // greedy algorithm to pick off the meta proteins with
  // most peptide mappings
  while (!peps_vector.empty()){
    sort(peps_vector.begin(), peps_vector.end(), setsAreEqualSize);
    pair<PeptideSet, MetaProtein> node = peps_vector.back();
    peps_vector.pop_back();
    if (node.first.size() == 0){ break; }// do not enter anything without peptide sizes
    result.insert(node);
    PeptideSet cur_peptides = node.first;
    // update the peptide sets for the rest of meta proteins
    for (vector< pair<PeptideSet, MetaProtein > >::iterator
           iter= peps_vector.begin();
         iter != peps_vector.end(); ++iter){
      PeptideSet peptides = (*iter).first;
      PeptideSet difference(Peptide::lessThan);
      set_difference(peptides.begin(), peptides.end(),
                     cur_peptides.begin(), cur_peptides.end(),
                     inserter(difference, difference.end()),
                    Peptide::lessThan);
      (*iter).first = difference;
    }
  }
  meta_mapping_ = result;
}

/**
 * Removes peptides from the map if the peptide
 * sequence belongs in more than one protein
 */
void SpectralCounts::makeUniqueMapping(){
  carp(CARP_INFO, "Filtering peptides that have more"
       "than one protein source");
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it){
    Peptide* peptide = it->first;
    int num_proteins = peptide->getNumPeptideSrc();
    if (num_proteins > 1){
      peptide_scores_.erase(it);
    }
  }

}


/**
 * \returns The name of the command.
 */
string SpectralCounts::getName() {
  return "spectral-counts";
}

/**
 * \returns The help string to be included in the crux usage statement.
 */
string SpectralCounts::getDescription() {
  return 
    "Rank proteins or peptides according to one of two quantification methods.";
}

COMMAND_T SpectralCounts::getCommand() {
  return SPECTRAL_COUNTS_COMMAND;
}

bool SpectralCounts::needsOutputDirectory() {
  return true;
}


// static comparison functions

/** 
 * Compare two sets of peptides and return true if the first unshared
 * peptide sequence in set one is lexically less than that in set
 * two.
 */
bool SpectralCounts::comparePeptideSets(PeptideSet set_one, 
                                        PeptideSet set_two){

  // compare each peptides in the two (sorted) sets
  PeptideSet::iterator iter1 = set_one.begin();
  PeptideSet::iterator iter2 = set_two.begin();

  while( iter1 != set_one.end() && iter2 != set_two.end() ){
    int diff = Peptide::triCompareSequence(*iter1, *iter2);
    if( diff < 0 ){
      return true;
    } else if(diff > 0 ){
      return false;
    } // else they are equal, compare the next

    ++iter1;
    ++iter2;
  }

  // all peptides were the same; are the sets the same size?
  if( set_one.size() == set_two.size() ){
    return false;
  } else if( set_one.size() > set_two.size() ){
    return false;
  } else { // one < two
    return true;
  }
}

/**
 * Comparison function for MetaProteins.  MetaProtein one is less than
 * MetaProtein two if the first non-matching protein id of one is less than
 * that of two.  
 * \returns True if one < two, false if one == two or one > two.
 */
bool SpectralCounts::compareMetaProteins(MetaProtein set_one, 
                                         MetaProtein set_two){

  // compare each protein in the two (sorted) sets
  MetaProtein::iterator iter1 = set_one.begin();
  MetaProtein::iterator iter2 = set_two.begin();

  while( iter1 != set_one.end() && iter2 != set_two.end() ){
    bool one_less_than_two = protein_id_less_than(*iter1, *iter2);
    bool two_less_than_one = protein_id_less_than(*iter1, *iter2);
    // different proteins one is less than the other
    if( one_less_than_two || two_less_than_one ){
      return one_less_than_two;
    }
    // else, they are the same, keep comparing
    ++iter1;
    ++iter2;
  }

  // all proteins were the same, are the sets the same size?
  if( set_one.size() == set_two.size() ){
    return false;
  } else if( set_one.size() > set_two.size() ){
    return false;
  } else { // one < two
    return true;
  }
}

/**
 * Compare the size of the PeptideSets in the two given pairs.
 * \returns True if the PeptideSets are the same size, else false.
 */
bool SpectralCounts::setsAreEqualSize(
  pair<PeptideSet, MetaProtein > peps_one ,
  pair<PeptideSet, MetaProtein > peps_two){
  return ((peps_one).first.size() < (peps_two).first.size());
}
