/**
 * \file PercolatorAdapter.cpp
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#include "PercolatorAdapter.h"

#include<map>
#include<string>

using namespace std;

/**
 * Constructor for PercolatorAdapter. 
 */
PercolatorAdapter::PercolatorAdapter() : Caller() {
  carp(CARP_DEBUG, "PercolatorAdapter::PercolatorAdapter");
  collection_ = new ProteinMatchCollection();
}

/**
 * Destructor for PercolatorAdapter
 */
PercolatorAdapter::~PercolatorAdapter() {
  carp(CARP_DEBUG, "PercolatorAdapter::~PercolatorAdapter");

  // delete match collections created by this adapter
  int collectionsDeleted = 0;
  for (vector<MatchCollection*>::iterator iter = match_collections_made_.begin();
       iter != match_collections_made_.end();
       ++iter) {
    delete *iter;
    ++collectionsDeleted;
  }
  // delete proteins created by this adapter
  int proteinsDeleted = 0;
  for (vector<PostProcessProtein*>::iterator iter = proteins_made_.begin();
       iter != proteins_made_.end();
       ++iter) {
    delete *iter;
    ++proteinsDeleted;
  }

  //delete collection_;  //TODO find out why this crashes
  carp(CARP_DEBUG, "PercolatorAdapter::~PercolatorAdapter - done. %d "
       "MatchCollections deleted.", collectionsDeleted);
}

/**
  * Calls Percolator's overridden Caller::writeXML_PSMs() and then
  * Collects all of the psm results
  */
void PercolatorAdapter::writeXML_PSMs() {
  carp(CARP_DEBUG, "PercolatorAdapter::writeXML_PSMs");
  Caller::writeXML_PSMs();
  addPsmScores();  
}

/**
 * Calls Percolator's overridden Caller::writeXML_Peptides() and then
 * Collects all of the peptide results from the fullset
 */
void PercolatorAdapter::writeXML_Peptides() {
  carp(CARP_DEBUG, "PercolatorAdapter::writeXML_Peptides");
  Caller::writeXML_Peptides();
  addPeptideScores();
}

/**
 * Calls Percolator's overriden Caller::writeXMLProteins() and then
 * Collects all of the protein results from fido
 */  
void PercolatorAdapter::writeXML_Proteins() {
  carp(CARP_DEBUG, "PercolatorAdapter::writeXML_Proteins");
  Caller::writeXML_Proteins();
  addProteinScores();
}

/**
 * Converts a set of Percolator scores into a Crux MatchCollection
 */
MatchCollection* PercolatorAdapter::psmScoresToMatchCollection() {

  // Create new MatchCollection object that will be the converted Percolator Scores
  MatchCollection* converted = new MatchCollection();
  match_collections_made_.push_back(converted);

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = fullset.begin();
    score_itr != fullset.end();
    score_itr++
    ) {

    bool is_decoy = score_itr->isDecoy();
    if (is_decoy) {
      continue;
    }

    PSMDescription* psm = score_itr->pPSM;
    int charge_state = parseChargeState(psm->id);
    Crux::Peptide* peptide = extractPeptide(psm, charge_state, is_decoy);

    SpectrumZState zState;
    zState.setSinglyChargedMass(psm->expMass, charge_state);
    // calcMass/expMass = singly charged mass
    Crux::Spectrum* spectrum = new Crux::Spectrum(
      psm->scan, psm->scan, zState.getMZ(), vector<int>(1, charge_state), ""
    );

    Crux::Match* match = new Crux::Match(peptide, spectrum, zState, is_decoy);
    match->setScore(PERCOLATOR_SCORE, score_itr->score);
    match->setScore(PERCOLATOR_QVALUE, psm->q);
    match->setScore(PERCOLATOR_PEP, psm->pep);

    converted->addMatch(match);
    match->setPostProcess(true); // so spectra get deleted when match does
    Crux::Match::freeMatch(match); // so match gets deleted when collection does
  }

  converted->forceScoredBy(PERCOLATOR_SCORE);
  converted->forceScoredBy(PERCOLATOR_QVALUE);
  converted->forceScoredBy(PERCOLATOR_PEP);
  converted->populateMatchRank(PERCOLATOR_SCORE);

  // sort by q-value
  converted->sort(PERCOLATOR_QVALUE);

  return converted;

}

/**
 * Adds PSM scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addPsmScores() {
  MatchCollection* matches = psmScoresToMatchCollection();
  collection_->addMatches(matches);
}

/**
 * Adds protein scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addProteinScores() {

  vector<ProteinMatch*> matches;
  map<const string,Protein> protein_scores = protEstimator->getProteins();
  
  for (
    map<const string,Protein>::iterator score_iter = protein_scores.begin();
    score_iter != protein_scores.end();
    score_iter++) {
    
    if (!score_iter->second.getIsDecoy()) {
      // Set scores
      ProteinMatch* protein_match = collection_->getProteinMatch(score_iter->second.getName());
      protein_match->setScore(PERCOLATOR_SCORE, -log(score_iter->second.getP()));
      protein_match->setScore(PERCOLATOR_QVALUE, score_iter->second.getQ());
      protein_match->setScore(PERCOLATOR_PEP, score_iter->second.getPEP());
      matches.push_back(protein_match);
    }
  }

  // set percolator score ranks
  std::sort(matches.begin(), matches.end(),
            PercolatorAdapter::comparePercolatorScores);
  int cur_rank = 0;
  for (vector<ProteinMatch*>::iterator iter = matches.begin();
       iter != matches.end();
       ++iter) {
    ProteinMatch* match = *iter;
    match->setRank(PERCOLATOR_SCORE, ++cur_rank);
  }

}

/**
 * Adds peptide scores from Percolator objects into a ProteinMatchCollection
 */
void PercolatorAdapter::addPeptideScores() {

  carp(CARP_DEBUG, "Setting peptide scores");

  // Iterate over each ScoreHolder in Scores object
  for (
    vector<ScoreHolder>::iterator score_itr = fullset.begin();
    score_itr != fullset.end();
    score_itr++
    ) {

    if (score_itr->isDecoy()) {
      continue;
    }

    PSMDescription* psm = score_itr->pPSM;
    string sequence;
    FLOAT_T peptide_mass;
    MODIFIED_AA_T* mod_seq = getModifiedAASequence(psm, sequence, peptide_mass);

    // Set scores
    PeptideMatch* peptide_match = collection_->getPeptideMatch(mod_seq);
    if (peptide_match == NULL) {
      carp(CARP_FATAL, "Cannot find peptide %s %i",
                       psm->getPeptideSequence().c_str(), score_itr->isDecoy());
    }
    peptide_match->setScore(PERCOLATOR_SCORE, score_itr->score);
    peptide_match->setScore(PERCOLATOR_QVALUE, psm->q);
    peptide_match->setScore(PERCOLATOR_PEP, psm->pep);

    free(mod_seq);

  }

}
  
/*
 *\returns the ProteinMatchCollection, to be called after Caller::run() is finished
 */
ProteinMatchCollection* PercolatorAdapter::getProteinMatchCollection() {
  return collection_;
}

/**
 * Given a Percolator psm_id in the form ".*_([0-9]+)_[^_]*",
 * find the charge state (matching group)
 */
int PercolatorAdapter::parseChargeState(
  string psm_id ///< psm to parse charge state from
) {
  size_t charge_begin, charge_end;

  charge_end = psm_id.rfind("_");
  if (charge_end < 0)
  {
    return -1;
  }

  charge_begin = psm_id.rfind("_", charge_end - 1) + 1;
  if (charge_begin < 0)
  {
    return -1;
  }

  string charge_str = psm_id.substr(charge_begin, charge_end - charge_begin);
  return atoi(charge_str.c_str());
}

/**
 * Compares two AbstractMatches by Percolator score
 */
bool PercolatorAdapter::comparePercolatorScores(
  AbstractMatch* lhs, ///< first match with Percolator scores to compare
  AbstractMatch* rhs ///< second match with Percolator scores to compare
) {
  if (!lhs->hasScore(PERCOLATOR_SCORE) || !rhs->hasScore(PERCOLATOR_SCORE)) {
    carp(CARP_FATAL, "Could not compare matches by Percolator score.");
  }
  return lhs->getScore(PERCOLATOR_SCORE) < rhs->getScore(PERCOLATOR_SCORE);
}

/**
* \returns a Crux peptide from the PSM
*/
Crux::Peptide* PercolatorAdapter::extractPeptide(
  PSMDescription* psm, ///< psm
  int charge_state, ///< charge state
  bool is_decoy ///< is psm a decoy?
  ) {

  string seq;
  FLOAT_T peptide_mass;
  
  MODIFIED_AA_T* mod_seq = getModifiedAASequence(psm, seq, peptide_mass);

  string full_peptide(psm->getFullPeptide());
  string n_term = "";
  string c_term = "";
  if (!full_peptide.empty()) {
    n_term += full_peptide[0];
    c_term += full_peptide[full_peptide.length() - 1];
  }

  PostProcessProtein* parent_protein = new PostProcessProtein();
  proteins_made_.push_back(parent_protein);
  parent_protein->setId((*psm->proteinIds.begin()).c_str());
  int start_idx = parent_protein->findStart(seq, n_term, c_term);

  Crux::Peptide* peptide = new Crux::Peptide(seq.length(), peptide_mass, parent_protein, start_idx);

  // add other proteins
  bool skip_one = true;
  for (set<string>::iterator iter = psm->proteinIds.begin();
       iter != psm->proteinIds.end();
       ++iter) {
    if (skip_one) {
      skip_one = false;
      continue;
    }
    PostProcessProtein* secondary_protein = new PostProcessProtein();
    proteins_made_.push_back(secondary_protein);
    secondary_protein->setId(iter->c_str());
    int secondary_idx = secondary_protein->findStart(seq, n_term, c_term);
    peptide->addPeptideSrc(
      new PeptideSrc(NON_SPECIFIC_DIGEST, secondary_protein, secondary_idx)
    );
  }

  peptide->setModifiedAASequence(mod_seq, is_decoy);

  free(mod_seq);
  return peptide;
}

/**
 * \returns the modified and unmodified peptide sequence
 * for the psm
 */
MODIFIED_AA_T* PercolatorAdapter::getModifiedAASequence(
  PSMDescription* psm, ///< psm -in
  string& seq, ///< sequence -out
  FLOAT_T& peptide_mass ///< calculated mass of peptide with modifications -out
  ) {

  std::stringstream ss_seq;
  string perc_seq = psm->getPeptideSequence();
  peptide_mass = 0.0;
  
  if (perc_seq.find("UNIMOD") != string::npos) {
    carp(CARP_FATAL, 
	 "UNIMOD modifications currently not supported:%s", 
	 perc_seq.c_str());
  }


  vector<pair<int, const AA_MOD_T*> > mod_locations_types;
  size_t count = 0;
  for (size_t seq_idx = 0; seq_idx < perc_seq.length(); seq_idx++) {
    if (perc_seq.at(seq_idx) == '[') {
      //modification found.
      size_t begin_idx = seq_idx+1;
      size_t end_idx = perc_seq.find(']', begin_idx);
      int mod_location = count-1;
      FLOAT_T delta_mass;

      from_string(delta_mass, perc_seq.substr(begin_idx, end_idx - begin_idx));
      carp(CARP_DEBUG,"seq:%s, i:%i m:%f", perc_seq.c_str(), seq_idx, delta_mass);
      peptide_mass += delta_mass;
      const AA_MOD_T* mod = get_aa_mod_from_mass(delta_mass);
      if (mod == NULL) {
	carp(CARP_FATAL, "Mod not found!");
      }

      mod_locations_types.push_back(pair<int, const AA_MOD_T*>(mod_location, mod));
      seq_idx = end_idx;
    } else {
      ss_seq << perc_seq.at(seq_idx);
      count++;
    }
  }
  seq = ss_seq.str();
  peptide_mass += Crux::Peptide::calcSequenceMass(seq.c_str(),
                  get_mass_type_parameter("isotopic-mass"));

  MODIFIED_AA_T* mod_seq;
  convert_to_mod_aa_seq(seq.c_str(), &mod_seq);
  for (size_t mod_idx = 0 ; mod_idx < mod_locations_types.size(); mod_idx++) {
    size_t seq_idx = mod_locations_types[mod_idx].first;
    const AA_MOD_T* mod = mod_locations_types[mod_idx].second;

    modify_aa(mod_seq+seq_idx, (AA_MOD_T*)mod);
  }
  return mod_seq;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

