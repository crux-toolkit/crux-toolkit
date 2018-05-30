#include <iterator>
#include "SpectralCounts.h"
#include "util/crux-utils.h"
#include "util/Params.h"
#include "io/OutputFiles.h"
#include "model/Peptide.h"
#include "model/ProteinPeptideIterator.h"
#include "io/SpectrumCollectionFactory.h"

using namespace std;
using namespace Crux;
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
    peptide_scores_unique_(Peptide::lessThan),
    peptide_scores_shared_(Peptide::lessThan),
    protein_scores_(protein_id_less_than),
    protein_scores_unique_(protein_id_less_than),
    protein_scores_shared_(protein_id_less_than),
    protein_supporting_peptides_(protein_id_less_than),
    protein_meta_protein_(protein_id_less_than),
    meta_mapping_(comparePeptideSets),
    meta_protein_scores_(compareMetaProteins),
    meta_protein_ranks_(compareMetaProteins) {
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
 * types of quantification: Raw counts (RAW), 
 * Normalized Spectral Abundance Factor (NSAF),
 * Distributed Normalized Spectral Abundance Factor (dNSAF),
 * EMPAI,
 * and Normalized Spectral Index (SIN). 
 * \returns 0 on successful completion.
 */
int SpectralCounts::main(int argc, char** argv) {
  getParameterValues(); // all the get_<type>_parameter calls here

  // open output files
  output_ = new OutputFiles(this);
  output_->writeHeaders();

  // get a set of matches that pass the threshold
  filterMatches();
  carp(CARP_INFO, "Number of matches passed the threshold %i", 
       matches_.size());

  if (matches_.empty()) {
    carp(CARP_FATAL, "No matches passed the threshold!");
  }

  // get a set of peptides
  getPeptideScores();
  if (unique_mapping_) {
    makeUniqueMapping();
  }
  carp(CARP_INFO, "Number of peptides %i", peptide_scores_.size());

  // quantify at either the peptide or protein level
  if (quantitation_ == PEPTIDE_QUANT_LEVEL) { // peptide level
    if (measure_ != MEASURE_RAW) {
      normalizePeptideScores();
    }
    writeRankedPeptides();

  } else if (quantitation_ == PROTEIN_QUANT_LEVEL) { // protein level
    
    getProteinScores();
    if (measure_ != MEASURE_RAW) {
      normalizeProteinScores();
      checkProteinNormalization();
    }
    carp(CARP_INFO, "Number of proteins %i", protein_scores_.size());
        
    if (parsimony_ != PARSIMONY_NONE) { //if parsimony is not none
      getProteinToPeptides();
      getMetaMapping();
      getProteinToMetaProtein();
      carp(CARP_INFO, "Number of meta proteins %i", meta_mapping_.size());

      if (parsimony_ == PARSIMONY_GREEDY) { //if parsimony is greedy
        performParsimonyAnalysis();
      }
      getMetaScores();
      getMetaRanks();
    }
    writeRankedProteins();
  } else {
    carp(CARP_FATAL, "Invalid quantification level.");
  }
  
  return 0;
}

/**
 * Collect all the parameter/option values we will need and store them
 * as member variables.
 */
void SpectralCounts::getParameterValues() {
  psm_file_ = Params::GetString("input PSMs");
  threshold_ = Params::GetDouble("threshold");
  database_name_ = Params::GetString("protein-database");
  unique_mapping_ = Params::GetBool("unique-mapping");
  quantitation_ = string_to_quant_level_type(Params::GetString("quant-level"));
  parsimony_ = string_to_parsimony_type(Params::GetString("parsimony"));
  measure_ = string_to_measure_type(Params::GetString("measure"));

  if (measure_ == MEASURE_SIN && Params::GetString("input-ms2").empty()) {
    carp(CARP_FATAL, "The SIN computation for spectral-counts requires "
                     "that the --input-ms2 option specify a file.");
  }

  bin_width_ = Params::GetDouble("mz-bin-width");
  
  threshold_type_ = get_threshold_type_parameter("threshold-type");
  custom_threshold_name_ = Params::GetString("custom-threshold-name");
  threshold_min_ = Params::GetBool("custom-threshold-min");
}

/**
 * For every protein that can be mapped from the set of 
 * peptides in PeptideToScore map, enter the protein and
 * the set of identified peptides it maps to, into 
 * ProteinToPeptide
 */
void SpectralCounts::getProteinToPeptides() {
  for (PeptideToScore::iterator pep_it = peptide_scores_.begin();
       pep_it != peptide_scores_.end(); ++pep_it) {
    Peptide* peptide = pep_it->first;
    for(PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
        iter!= peptide->getPeptideSrcEnd();
        ++iter) {
      PeptideSrc* peptide_src = *iter; 
      Protein* protein = peptide_src->getParentProtein();
      if (protein_supporting_peptides_.find(protein) == 
          protein_supporting_peptides_.end()) {
        PeptideSet newset(Peptide::lessThan);
        protein_supporting_peptides_.insert(make_pair(protein, newset));
      }
      protein_supporting_peptides_[protein].insert(peptide);
    }
  }
}


/**
 * Enters the mapping of protein to its metaProtein
 * into ProteinToMetaProtein. MetaProteins are retreieved
 * from MetaMapping
 *
 */
void SpectralCounts::getProteinToMetaProtein() {
  // for every meta protein
  for (MetaMapping::iterator meta_protein_it = meta_mapping_.begin();
       meta_protein_it != meta_mapping_.end(); ++meta_protein_it) {
    MetaProtein proteins = meta_protein_it->second;
    // for every protein in the meta protein
    for (MetaProtein::iterator proteins_it = proteins.begin();
         proteins_it != proteins.end(); ++proteins_it) {
      // create a mapping of protein to meta protein
      protein_meta_protein_.insert(make_pair((*proteins_it), proteins));
    }
  }
}

/**
 * calculates the protein scores based upon the dNSAF metric.  First,
 * the spectral counts for the peptides unique to each protein is added,
 * then the spectral counts for the peptide shared amongst the proteins
 * are added with a distribution factor based upon the number of unique
 * counts that each protein is assigned
 */
void SpectralCounts::getProteinScoresDNSAF() {

  //calculate unique scores (unique peptides).
  for (PeptideToScore::iterator pep_it = peptide_scores_unique_.begin();
    pep_it != peptide_scores_unique_.end();
    ++pep_it) {

    Peptide* peptide = pep_it->first;
    FLOAT_T pep_score = pep_it->second;
    for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
         iter != peptide->getPeptideSrcEnd();
         ++iter) {
      PeptideSrc* peptide_src = *iter;
      Protein* protein = peptide_src->getParentProtein();
      if (protein_scores_unique_.find(protein) == protein_scores_unique_.end()) {
        protein_scores_unique_[protein] = 0;
        protein_scores_shared_[protein] = 0;
      }

      protein_scores_unique_[protein] += pep_score;
    }   
  }

  //Handle shared peptides
  //For each peptide that is shared across multiple proteins,
  //we calculate a distribution factor based upon the unique
  //peptides assigned to each protein that is sharing this peptide
  //we then add the shared peptide count * the protein's distribution
  //to each protein that is sharing that peptide.
  //There is a wierd feature that a protein with no unique peptide will
  //have a dNSAF score of zero... I don't know if we should include the
  //zero score in the list of proteins...
  for (PeptideToScore::iterator pep_it = peptide_scores_shared_.begin();
    pep_it != peptide_scores_shared_.end();
    ++pep_it) {
    Peptide* peptide = pep_it->first;
    FLOAT_T shared_pep_score = pep_it->second;

    double unique_sum = 0.0;
   
    for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
         iter != peptide->getPeptideSrcEnd();
         ++iter) {
      PeptideSrc* peptide_src = *iter;
      
      Protein* protein = peptide_src->getParentProtein();
      if (protein_scores_unique_.find(protein) != protein_scores_unique_.end()) {
        unique_sum += protein_scores_unique_[protein];
      }
        
    }   
  
    if (unique_sum != 0) {
      for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
           iter != peptide->getPeptideSrcEnd();
           ++iter) {
        PeptideSrc* peptide_src = *iter;
        Protein* protein = peptide_src->getParentProtein();
        if (protein_scores_unique_.find(protein) != protein_scores_unique_.end()) {
          FLOAT_T d_factor = protein_scores_unique_[protein] / unique_sum;
          protein_scores_shared_[protein] += d_factor * shared_pep_score;
        }
      }
    }
  }

  //add up shared and unique scores to get final score
  for (ProteinToScore::iterator prot_iter = protein_scores_unique_.begin();
    prot_iter != protein_scores_unique_.end();
    ++prot_iter) {
    Protein* protein = prot_iter->first;
    FLOAT_T unique_score = prot_iter->second;
    protein_scores_[protein] = unique_score;
  }

  for (ProteinToScore::iterator prot_iter = protein_scores_shared_.begin();
    prot_iter != protein_scores_shared_.end();
    ++prot_iter) {
    Protein *protein = prot_iter->first;
    FLOAT_T shared_score = prot_iter->second;
    protein_scores_[protein] += shared_score;
  }
}

/**
 * A score for each protein is calculated by summing 
 * the scores of each peptide that belongs to a protein
 */
void SpectralCounts::getProteinScores() {

  if (measure_ == MEASURE_DNSAF) {
    getProteinScoresDNSAF();
  } else {

    // iterate through each peptide
    for (PeptideToScore::iterator pep_it = peptide_scores_.begin();
         pep_it != peptide_scores_.end(); ++pep_it) {
      Peptide* peptide = pep_it->first;
      FLOAT_T pep_score = pep_it->second;
      for (PeptideSrcIterator iter = peptide->getPeptideSrcBegin();
           iter != peptide->getPeptideSrcEnd();
           ++iter) {
        PeptideSrc* peptide_src =*iter;
        Protein* protein = peptide_src->getParentProtein();
        if (protein_scores_.find(protein) == protein_scores_.end()) {
          protein_scores_.insert(make_pair(protein, 0.0));
        }
        protein_scores_[protein] += pep_score;
      }
    }
  }
}

/**
 * Takes the PeptideToScores map and updates all
 * values with normalized values. Normalized by sum
 * of all scores and then by the peptide length
 *
 */
void SpectralCounts::normalizePeptideScores() {
  carp(CARP_DEBUG, "Normalizing peptide scores");
  FLOAT_T total = 0.0;

  // calculate sum of all scores
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it) {
    FLOAT_T score = it->second;
    total += score;
  }

  // normalize by sum of scores and length
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it) {
    FLOAT_T score = it->second;
    Peptide* peptide = it->first;
    it->second = score / total / (FLOAT_T)peptide->getLength();

  }

}


/**
 * Changes the scores in protien_scores_ to either be divided by the
 * sum of all scores times the peptide length (SIN, NSAF) or to be the
 * final emPAI score.
 */
void SpectralCounts::normalizeProteinScores() {
  if (measure_ == MEASURE_EMPAI) {
    computeEmpai();
  } 

  carp(CARP_DEBUG, "Normalizing protein scores");
  FLOAT_T total = 0.0;
    
  // calculate sum of all scores
  for (ProteinToScore::iterator it = protein_scores_.begin();
       it != protein_scores_.end(); ++it) {

    FLOAT_T score = it->second;
    Protein* protein = it->first;
    if ( (measure_ == MEASURE_NSAF) || (measure_ == MEASURE_DNSAF)) {
      score = score / (FLOAT_T)protein->getLength();
    }
    total += score;
  }
    
  // normalize by sum of all scores
  for (ProteinToScore::iterator it = protein_scores_.begin();
       it != protein_scores_.end(); ++it) {

    FLOAT_T score = it->second;
    Protein* protein = it->first;
    score = score / total;
    // normalize by length
    if (measure_ == MEASURE_NSAF || 
        measure_ == MEASURE_DNSAF || 
        measure_ == MEASURE_SIN ) {
      score = score / (FLOAT_T)protein->getLength();
    }
    it->second = score;
  }
}

/**
 * Checks that the normalized scores add up to one
 */
void SpectralCounts::checkProteinNormalization() {

  FLOAT_T sum = 0;

  for (ProteinToScore::iterator iter = protein_scores_.begin();
    iter != protein_scores_.end();
    ++iter) {

    FLOAT_T score = iter->second;
    if (measure_ == MEASURE_SIN) {
      //The normalized values of sin do not add up to 1, but they
      //should if you multiply the length back in...
      Protein* protein = iter->first;
      score = score * (FLOAT_T)protein->getLength();
    }

    sum += score;
  }

  if (fabs(sum-1.0) > 0.0001) {
    carp(CARP_WARNING, "Normalized protein scores do not add up to one!:%f", sum);
    for (ProteinToScore::iterator iter = protein_scores_.begin();
      iter != protein_scores_.end();
      ++iter) {

      carp(CARP_DEBUG, "%s %f", iter->first->getIdPointer().c_str(), iter->second);
    }
  }
}

/**
 * Computes the 10^(observed/total) - 1 score for each protein.  Assumes
 * that protein_scores_ has been populated with the count of observed
 * unique peptides for each protein.
 */
void SpectralCounts::computeEmpai() {
  PeptideConstraint* constraint = PeptideConstraint::newFromParameters();

  ProteinToScore::iterator it = protein_scores_.begin();;
  for (; it != protein_scores_.end(); ++it) {
    FLOAT_T observed_peptides = it->second;
    Protein* protein = it->first;
    ProteinPeptideIterator* iter = new ProteinPeptideIterator(protein, 
                                                              constraint);
    FLOAT_T possible_peptides = iter->getTotalPeptides();

    it->second = pow(10, (observed_peptides / possible_peptides)) - 1.0;

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
                                        Crux::SpectrumCollection* spectra) {
  FLOAT_T match_intensity = 0;
  char* peptide_seq = match->getSequence();
  MODIFIED_AA_T* modified_sequence = match->getModSequence();
  int charge = match->getCharge();
  Spectrum* temp = match->getSpectrum();
  int scan = temp->getFirstScan();

  Spectrum* spectrum = spectra->getSpectrum(scan);

  if (spectrum == NULL) {
    carp(CARP_FATAL, "scan: %d doesn't exist or not found!");
    return 0.0;
  }

  Ion* ion;
  SCORER_TYPE_T score_type = XCORR;
  IonConstraint* ion_constraint =
    IonConstraint::newIonConstraintSmart(score_type, charge);
  IonSeries* ion_series = new IonSeries(ion_constraint, charge);
  ion_series->update(peptide_seq, modified_sequence);
  ion_series->predictIons();
  for (IonIterator ion_it = ion_series->begin();
       ion_it != ion_series->end(); ++ion_it) {
    ion = (*ion_it);
    if (ion -> getType() == B_ION || ion -> getType() == Y_ION) {
      if (!ion->isModified()) {
        Peak * peak = spectrum->getNearestPeak(ion->getMassZ(),
                                                bin_width_);
        if (peak != NULL) {
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
void SpectralCounts::getPeptideScores() {
  Crux::SpectrumCollection* spectra = NULL;

  // for SIN, parse out spectrum collection from ms2 fiel
  if( measure_ == MEASURE_SIN ) {
    spectra = SpectrumCollectionFactory::create(Params::GetString("input-ms2"));
  }

  for(set<Match*>::iterator match_it = matches_.begin();
      match_it != matches_.end(); ++match_it) {

    FLOAT_T match_intensity = 1; // for NSAF just count each for the peptide/

    Match* match = (*match_it);
    // for sin, calculate total ion intensity for match by
    // summing up peak intensities
    if (measure_ == MEASURE_SIN) {
      match_intensity = sumMatchIntensity(match, spectra);
    }

    // add ion_intensity to peptide scores
    Peptide* peptide = match->getPeptide();

    if ( measure_ == MEASURE_DNSAF) {
      //keep track of unique and shared mappings.
      if (peptide->getNumPeptideSrc() > 1) {
        if (peptide_scores_shared_.find(peptide) == peptide_scores_shared_.end()) {
          peptide_scores_shared_.insert(make_pair(peptide, 0.0));
        }
        peptide_scores_shared_[peptide] += match_intensity;
      } else {
        if (peptide_scores_unique_.find(peptide) == peptide_scores_unique_.end()) {
          peptide_scores_unique_.insert(make_pair(peptide, 0.0));
        }
        peptide_scores_unique_[peptide] += match_intensity;
      }
    }

    if (peptide_scores_.find(peptide) ==  peptide_scores_.end()) {
      peptide_scores_.insert(make_pair(peptide, 0.0));
    }
    peptide_scores_[peptide] += match_intensity;

  }

  if (measure_ == MEASURE_SIN) {
    delete spectra;
  }

  // for emPAI we just need a count of unique peptides
  if (measure_ == MEASURE_EMPAI) {
    PeptideToScore::iterator itr = peptide_scores_.begin();
    for(; itr != peptide_scores_.end(); ++itr) {
      itr->second = 1; // count each unique peptide once
    }
  }
}


/**
 * Create a set of matches, all with an XCORR rank == 1 and all of which
 * have a qvalue score lower than user-specified threshold.
 */
void SpectralCounts::filterMatches() {

  match_collection_ = parser_.create(
    psm_file_.c_str(),
    database_name_.c_str());
  carp(CARP_INFO, "Number of matches:%d", match_collection_->getMatchTotal());

  switch(threshold_type_) {
    case THRESHOLD_NONE:
      filterMatchesNone();
      break;
    case THRESHOLD_QVALUE:
      filterMatchesQValue();
      break;
    case THRESHOLD_CUSTOM:
      filterMatchesCustom();
      break;
    case THRESHOLD_INVALID:
    case NUMBER_THRESHOLD_TYPES:
      carp(CARP_FATAL, "Invalid threshold type");
  }

}

void SpectralCounts::filterMatchesNone() {

  MatchIterator match_iterator(match_collection_);

  while(match_iterator.hasNext()) {
    Match* current_match = match_iterator.next();
    if (!current_match->isDecoy()) {
      matches_.insert(current_match);
    }
  }

}

void SpectralCounts::invalidCustomScore() {

    ostringstream oss;

    oss << "Need valid name for custom threshold.  Available ones are:";

    //print out custom names
    for (int idx=0;idx < NUMBER_SCORER_TYPES;idx++) {
    
      if (match_collection_->getScoredType((SCORER_TYPE_T)idx)) {
        oss << endl << "\t" << scorer_type_to_string((SCORER_TYPE_T)idx);
      }

    }
    

    vector<string> custom_score_names;
    match_collection_->getCustomScoreNames(custom_score_names);

    for (int idx=0;idx < custom_score_names.size();idx++) {
      oss << endl << "\t" << custom_score_names[idx];
    }

    string die_str = oss.str();
    carp(CARP_FATAL, die_str);
}


/**
 * filters matches based upon a custom threshold (not q-value)
 */
void SpectralCounts::filterMatchesCustom() {

  if (custom_threshold_name_.empty()) {
    invalidCustomScore();
  }


  //first try to map from tab delimited headers
  SCORER_TYPE_T scorer;
  bool tab_header = string_to_scorer_type(custom_threshold_name_.c_str(), &scorer);

  if (tab_header) {
    //must be a custom field, this can happen with pep xml or mzid.
    filterMatchesScore(scorer);    
  } else {
    filterMatchesCustomScore();
  }


}

/**
 * filters matches based upon a SCORER_TYPE_T
 */
void SpectralCounts::filterMatchesScore(
  SCORER_TYPE_T scorer ///< scorer to use
  ) {

  MatchIterator match_iterator(match_collection_, scorer, false);
  
  while (match_iterator.hasNext()) {
    Match* match = match_iterator.next();
    if (!match->isDecoy()) {
      if (threshold_min_) {
        if (match->getScore(scorer) <= threshold_) {
          matches_.insert(match);
        } 
      } else {
        if (match->getScore(scorer) >= threshold_) {
          matches_.insert(match);
        }
      }
    }
  }
}

/**
 * filters matches based upon a custom score that is not SCORER_TYPE_T
 */
void SpectralCounts::filterMatchesCustomScore() {

  MatchIterator match_iterator(match_collection_);

  while(match_iterator.hasNext()) {
    Match* match = match_iterator.next();
    if (!match->isDecoy()) {
      FLOAT_T score;
      bool success = match->getCustomScore(custom_threshold_name_, score);
  
      if (!success) {
        invalidCustomScore();
      }
      if (threshold_min_) {
        if (score <= threshold_) {
          matches_.insert(match);
        }
      } else {
        if (score >= threshold_) {
          matches_.insert(match);
        }
      }
    }
  }
}

void SpectralCounts::filterMatchesQValue() {
  //assume we are using 
  // figure out which qvalue we are using
  SCORER_TYPE_T qval_type = get_qval_type(match_collection_);
  if (qval_type == INVALID_SCORER_TYPE) {
    carp(CARP_FATAL, "The matches in %s do not have q-values.\n"
                     "Please either provide a file with q-values or "
                     "set threshold-type to \"none\" or \"custom.\"",
    psm_file_.c_str());
  }

  carp(CARP_DETAILED_DEBUG,
    "filterMatches(): Getting match iterator for %s", 
    scorer_type_to_string(qval_type));

  MatchIterator* match_iterator = new MatchIterator(match_collection_, qval_type, true);

  while (match_iterator->hasNext()) {
    Match* match = match_iterator->next();
     
    carp(CARP_DEBUG, "xcorr rank:%d q-value:%f", match->getRank(XCORR), match->getScore(qval_type));
    if (match->isDecoy()) {
      continue;
    } 
    if ((qval_type == DECOY_XCORR_QVALUE || qval_type == LOGP_QVALUE_WEIBULL_XCORR) && 
      (match->getRank(XCORR) != 1)) {
      continue;
    }
      
    // find a qvalue score lower than threshold
    if (match->getScore(qval_type) != FLT_MIN &&
        match->getScore(qval_type) <= threshold_)  {
      matches_.insert(match);
    }
  } // next match
  delete match_iterator;
}

/**
 * Figures out which kind of q-value was scored for this match collection.
 * \returns PERCOLATOR_QVALUE, QRANKER_QVALUE, or DECOY_XCORR_QVALUE
 * if any of those were scored or INVALID_SCORER_TYPE if none were scored. 
 */
SCORER_TYPE_T SpectralCounts::get_qval_type(
  MatchCollection* match_collection) {
  SCORER_TYPE_T scored_type = INVALID_SCORER_TYPE;

  if (match_collection->getScoredType(PERCOLATOR_QVALUE)) {
    scored_type =  PERCOLATOR_QVALUE;
  } else if (match_collection->getScoredType(QRANKER_QVALUE)) {
    scored_type = QRANKER_QVALUE;
  } else if (match_collection->getScoredType(DECOY_XCORR_QVALUE)) {
    scored_type = DECOY_XCORR_QVALUE;
  } else if(match_collection->getScoredType(BARISTA_QVALUE)) {
    scored_type = BARISTA_QVALUE;
  }

  return scored_type;
}

void SpectralCounts::writeRankedPeptides() {
  // rearrange pairs to sort by score
  vector<pair<FLOAT_T, Peptide*> > scoreToPeptide;
  for (PeptideToScore::iterator it = peptide_scores_.begin();
       it != peptide_scores_.end(); ++it) {
    scoreToPeptide.push_back(make_pair(it->second, it->first));
  }
  
  sort(scoreToPeptide.begin(), scoreToPeptide.end(), sortRankedPeptides);
  output_->writeRankedPeptides(scoreToPeptide);
}

void SpectralCounts::writeRankedProteins() {
  bool isParsimony = !protein_meta_protein_.empty();
  // reorganize the protein,score pairs to sort by score
  vector<boost::tuple<FLOAT_T, Protein*, int> > proteins;
  for (ProteinToScore::iterator it = protein_scores_.begin(); 
       it != protein_scores_.end(); ++it) {
    int rank = -1;
    if (isParsimony) {
      /* This find doesn't seem to work for all proteins within a metaprotein, 
         needs to be debugged
         if we are going to use it. For now, do a brute force find. (SJM _2018_05_26)   
      MetaToRank::const_iterator lookup =
        meta_protein_ranks_.find(protein_meta_protein_[it->first]);
      if (lookup != meta_protein_ranks_.end()) {
        rank = lookup->second;
      } else {
      */
      for (MetaToRank::iterator iter = meta_protein_ranks_.begin();
        iter != meta_protein_ranks_.end();
	++iter) {
	  MetaProtein proteins = iter->first;
	  for (MetaProtein::iterator protein_it = proteins.begin();
	       protein_it != proteins.end(); ++protein_it) {
	    Protein* protein = (*protein_it);
	    if (protein->getId() == it->first->getId()) {
	      carp(CARP_DEBUG, "Found protein %s",protein->getId().c_str());
	      rank = iter->second;
	    }
	  }
      }
      // }
    }
    proteins.push_back(boost::make_tuple(it->second, it->first, rank));
  }

  sort(proteins.begin(), proteins.end(), sortRankedProteins);
  output_->writeRankedProteins(proteins, isParsimony);
}

bool SpectralCounts::sortRankedPeptides(
  const pair<FLOAT_T, Peptide*>& x,
  const pair<FLOAT_T, Peptide*>& y) {
  return x.first != y.first ? x.first > y.first : !Peptide::lessThan(y.second, x.second);
}

bool SpectralCounts::sortRankedProteins(
  const boost::tuple<FLOAT_T, Protein*, int>& x,
  const boost::tuple<FLOAT_T, Protein*, int>& y) {
  FLOAT_T xScore = x.get<0>();
  FLOAT_T yScore = y.get<0>();
  if (xScore != yScore) {
    return xScore > yScore;
  }
  int xRank = x.get<2>();
  int yRank = y.get<2>();
  if (xRank != yRank) {
    return xRank <= yRank;
  }
  return y.get<1>()->getIdPointer().compare(x.get<1>()->getIdPointer()) > 0;
}

/**
 * Fills in the MetaMapping with entries of set of 
 * peptides that can be found in every protein in
 * the meta protein
 *
 */
void SpectralCounts::getMetaMapping() {
  carp(CARP_DEBUG, "Creating a mapping of meta protein to peptides");
  int count = 0;
  for(ProteinToPeptides::iterator prot_it = protein_supporting_peptides_.begin();
       prot_it != protein_supporting_peptides_.end(); ++prot_it) {
    Protein* protein = prot_it->first;
    PeptideSet pep_set = prot_it->second;

    if (meta_mapping_.find(pep_set) == meta_mapping_.end()) {
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
void SpectralCounts::getMetaScores() {
  carp(CARP_DEBUG, "Finding scores of meta proteins");
  for (MetaMapping::iterator meta_it = meta_mapping_.begin();
       meta_it != meta_mapping_.end(); ++meta_it) {
    MetaProtein proteins = (*meta_it).second;
    FLOAT_T top_score = -1.0;
    carp(CARP_DEBUG, "Meta protein");
    for (MetaProtein::iterator protein_it = proteins.begin();
         protein_it != proteins.end(); ++protein_it) {
      Protein* protein = (*protein_it);
      FLOAT_T score = protein_scores_[protein];
      top_score = max(score, top_score);
      carp(CARP_DEBUG, "   Protein %s score:%g", protein->getIdPointer().c_str(), score);
    }
    meta_protein_scores_.insert(make_pair(proteins, top_score));
  }

}

/**
 * Takes a mapping of MetaProteins to scores and returns 
 * a mapping of set of peptides to rank
 *
 */
void SpectralCounts::getMetaRanks() {
  carp(CARP_DEBUG, "Finding ranks of meta proteins");
  vector< pair<FLOAT_T, MetaProtein> > metaVector;
  for (MetaToScore::iterator meta_it = meta_protein_scores_.begin();
       meta_it != meta_protein_scores_.end(); ++meta_it) {
    MetaProtein proteins = (*meta_it).first;
    FLOAT_T score = (*meta_it).second;
    metaVector.push_back(make_pair(score, proteins));
  }
  sort(metaVector.begin(), metaVector.end(), compareMetaScorePair);

  int cur_rank = 1;
  FLOAT_T last_score = -1;
  for (size_t idx = 0;idx < metaVector.size();idx++) {
    MetaProtein proteins = metaVector[idx].second;
    FLOAT_T cur_score = metaVector[idx].first;
    if (cur_score != last_score) {
      cur_rank = idx+1;
    }

    carp(CARP_DEBUG, "Meta Protein score:%g rank:%i",cur_score,cur_rank);
     for (MetaProtein::iterator protein_it = proteins.begin();
         protein_it != proteins.end(); ++protein_it) {
      Protein* protein = (*protein_it);
      FLOAT_T score = protein_scores_[protein];
      carp(CARP_DEBUG, "   Protein %s score:%g", protein->getIdPointer().c_str(), score);
    }
    
    meta_protein_ranks_.insert(make_pair(proteins, cur_rank));
    last_score = cur_score;
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
void SpectralCounts::performParsimonyAnalysis() {
  carp(CARP_DEBUG, "Performing Greedy Parsimony analysis");
  MetaMapping result(comparePeptideSets);
  vector< pair<PeptideSet, MetaProtein > > peps_vector;

  // get all meta mappings into a vector 
  for (MetaMapping::iterator meta_iter = meta_mapping_.begin();
       meta_iter != meta_mapping_.end(); ++meta_iter) {
    peps_vector.push_back((*meta_iter));
  }


  // greedy algorithm to pick off the meta proteins with
  // most peptide mappings
  while (!peps_vector.empty()) {
    sort(peps_vector.begin(), peps_vector.end(), setsAreEqualSize);
    pair<PeptideSet, MetaProtein> node = peps_vector.back();
    peps_vector.pop_back();
    if (node.first.size() == 0) { break; }// do not enter anything without peptide sizes
    result.insert(node);
    PeptideSet cur_peptides = node.first;
    // update the peptide sets for the rest of meta proteins
    for (vector< pair<PeptideSet, MetaProtein > >::iterator
         iter = peps_vector.begin();
         iter != peps_vector.end(); ++iter) {
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
void SpectralCounts::makeUniqueMapping() {

  if (measure_ == MEASURE_DNSAF) {
    carp(CARP_WARNING, "--unique-mapping ignored for dNSAF!");
  } else {
  carp(CARP_DEBUG, "Filtering peptides that have more"
         "than one protein source");
    for (PeptideToScore::iterator it = peptide_scores_.begin();
         it != peptide_scores_.end(); ++it) {
      Peptide* peptide = it->first;
      int num_proteins = peptide->getNumPeptideSrc();
      if (num_proteins > 1) {
        peptide_scores_.erase(it);
      }
    }
  }
}


/**
 * \returns The name of the command.
 */
string SpectralCounts::getName() const {
  return "spectral-counts";
}

/**
 * \returns The help string to be included in the crux usage statement.
 */
string SpectralCounts::getDescription() const {
  return
    "[[nohtml:Quantify peptides or proteins using one of three spectral "
    "counting methods.]]"
    "[[html:<p>Given a collection of scored PSMs, produce a list of proteins "
    "or peptides ranked by a quantification score. Spectral-counts supports "
    "four types of quantification: Normalized Spectral Abundance Factor "
    "(NSAF), Distributed Normalized Spectral Abundance (dNSAF), Normalized "
    "Spectral Index (SI<sub>N</sub>) and Exponentially Modified Protein "
    "Abundance Index (emPAI). The NSAF method is from <a href=\""
    "http://www.ncbi.nlm.nih.gov/pubmed/17138671\">Paoletti et "
    "al. (2006)</a>. The SI<sub>N</sub> method is from <a href=\""
    "http://www.nature.com/nbt/journal/v28/n1/abs/nbt.1592.html\">Griffin et "
    "al. (2010)</a>. The emPAI method was first described in <a href=\""
    "http://www.mcponline.org/content/4/9/1265\">Ishihama et al (2005)</a>. "
    "The quantification methods are defined below and in the following "
    "paper:<blockquote>S McIlwain, M Mathews, M Bereman, EW Rubel, MJ "
    "MacCoss, and WS Noble.  "
    "<a href=\"http://www.biomedcentral.com/1471-2105/13/308/abstract\">"
    "\"Estimating relative abundances of proteins from shotgun proteomics "
    "data.\"</a>  <em>BMC Bioinformatics</em>. 13:308, 2012.</blockquote>"  
    "</p><h3>Protein Quantification</h3>"
    "<ol><li>For each protein in a given database, the NSAF "
    "score is:<br>$$NSAF_N=\\frac{S_N/L_N}{\\sum_{i=1}^ns_i/L_i}$$<br>where:"
    "<ul><li>N is protein index</li><li>S<sub>N</sub> is the number of peptide "
    "spectra matched to the protein</li><li>L<sub>N</sub> is the length of "
    "protein N</li><li>n is the total number of proteins in the input "
    "database</li></ul></li><li>For each protein in a given database, the dNSAF "
    "score is:<br>$$NSAF_N=\\frac{\\frac{uSpc_N+(d)sSpc_N}{uL_N+sL_N}}{\\frac"
    "{uSpc_i+(d)sSpc_i}{uL_i+sL_i}}$$<br>where:<ul><li>N is the protein "
    "index</li><li>uSpc<sub>N</sub> is the unique number spectra matched to the "
    "protein index</li><li>sSpc<sub>N</sub> is the shared number peptide "
    "spectra matched to the protein index</li><li>L<sub>N</sub> is the length "
    "of protein N</li><li>n is the total number of proteins in the input "
    "database</li><li>d is the distribution factor of peptide K to protein N, "
    "given by<br>$$d=\\frac{uSpc_N}{\\sum_{i=1}^nuSpc_i}$$</li></ul></li><li>"
    "For each protein in a given database, the SI<sub>N</sub> score is:<br>"
    "$$SI_N=\\frac{\\sum_{j=1}^{p_N}(\\sum_{k=1}^{s_j}i_k)}{L_N(\\sum_{j=1}^n"
    "SI_j)}$$<br>where:<ul><li>N is protein index</li><li>p<sub>n</sub> is the "
    "number of unique peptides in protein N</li><li>s<sub>j</sub> is the "
    "number of spectra assigned to peptide j</li><li>i<sub>k</sub> is the "
    "total fragment ion intensity of spectrum k</li><li>L<sub>N</sub> is the "
    "length of protein N</li></ul></li><li>For each protein in a given database, "
    "the emPAI score is:<br>$$emPAI=10^{\\frac{N_{observed}}{N_{observable}}}-1"
    "$$<br>where:<ul><li>N<sub>observed</sub> is the number of experimentally "
    "observed peptides with scores above a specified threshold.</li><li>N<sub>"
    "observable</sub> is the calculated number of observable peptides for the "
    "protein given the search constraints.</li></ul></li></ol><h3>Peptide "
    "Quantification</h3><ol><li>For each peptide in a given database, the NSAF "
    "score is:<br>$$NSAF_N=\\frac{S_N/L_N}{\\sum_{i=1}^ns_i/L_i}$$<br>where: "
    "<ul><li>N is the peptide index</li><li>S<sub>N</sub> is the number "
    "spectra matched to peptide N</li><li>L<sub>N</sub> is the length of "
    "peptide N</li><li>n is the total number of peptides in the input "
    "database</li></ul></li><li>For each peptide in a given database, the "
    "SI<sub>N</sub> score is:<br>$$SI_N=\\frac{(\\sum_{k=1}^{S_N}i_k)}{L_N("
    "\\sum_{j=1}^nSI_J)}$$<br>where:<ul><li>N is the peptide index</li><li>"
    "S<sub>N</sub> is the number of spectra assigned to peptide N</li><li>"
    "i<sub>k</sub> is the total fragment ion intensity of spectrum k</li><li>"
    "L<sub>N</sub> is the length of peptide N</li></ul></li></ol>]]";
}

vector<string> SpectralCounts::getArgs() const {
  string arr[] = {
    "input PSMs"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector<string> SpectralCounts::getOptions() const {
  string arr[] = {
    "verbosity",
    "parameter-file",
    "parsimony",
    "threshold",
    "threshold-type",
    "input-ms2",
    "spectrum-parser",
    "fileroot",
    "output-dir",
    "overwrite",
    "unique-mapping",
    "quant-level",
    "measure",
    "custom-threshold-name",
    "custom-threshold-min",
    "mzid-use-pass-threshold",
    "protein-database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

vector< pair<string, string> > SpectralCounts::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("spectral-counts.target.txt",
    "a tab-delimited text file containing the protein IDs and their "
    "corresponding scores, in sorted order."));
  outputs.push_back(make_pair("spectral-counts.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair("spectral-counts.log.txt",
    "All messages written to standard error."));
  return outputs;
}

COMMAND_T SpectralCounts::getCommand() const {
  return SPECTRAL_COUNTS_COMMAND;
}

bool SpectralCounts::needsOutputDirectory() const {
  return true;
}


// static comparison functions

/** 
 * Compare two sets of peptides and return true if the first unshared
 * peptide sequence in set one is lexically less than that in set
 * two.
 */
bool SpectralCounts::comparePeptideSets(PeptideSet set_one, 
                                        PeptideSet set_two) {

  // compare each peptides in the two (sorted) sets
  PeptideSet::iterator iter1 = set_one.begin();
  PeptideSet::iterator iter2 = set_two.begin();

  while( iter1 != set_one.end() && iter2 != set_two.end() ) {
    int diff = Peptide::triCompareSequence(*iter1, *iter2);
    if (diff != 0) {
      return diff < 0;
    }
    // else they are equal, compare the next
    ++iter1;
    ++iter2;
  }

  // all peptides were the same; are the sets the same size?
  return set_one.size() < set_two.size();
}

/**
 * Comparison function for MetaProteins.  MetaProtein one is less than
 * MetaProtein two if the first non-matching protein id of one is less than
 * that of two.  
 * \returns True if one < two, false if one == two or one > two.
 */
bool SpectralCounts::compareMetaProteins(MetaProtein set_one, 
                                         MetaProtein set_two) {
  // compare each protein in the two (sorted) sets
  MetaProtein::iterator iter1 = set_one.begin();
  MetaProtein::iterator iter2 = set_two.begin();

  while (iter1 != set_one.end() && iter2 != set_two.end()) {
    // different proteins one is less than the other
    if (protein_id_less_than(*iter1, *iter2)) {
      return true;
    }
    // else, they are the same, keep comparing
    ++iter1;
    ++iter2;
  }

  // all proteins were the same, are the sets the same size?
  return set_one.size() < set_two.size();
}

/**
 * Compare the size of the PeptideSets in the two given pairs.
 * \returns True if the PeptideSets are the same size, else false.
 */
bool SpectralCounts::setsAreEqualSize(
  const pair<PeptideSet, MetaProtein>& peps_one ,
  const pair<PeptideSet, MetaProtein>& peps_two) {
  return ((peps_one).first.size() < (peps_two).first.size());
}

bool SpectralCounts::compareMetaScorePair(
  const std::pair<FLOAT_T, MetaProtein>& x,
  const std::pair<FLOAT_T, MetaProtein>& y) {
  return (x.first != y.first) ? x.first > y.first : !compareMetaProteins(x.second, y.second);
}

