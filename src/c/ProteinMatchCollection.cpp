/**
 * \file ProteinMatchCollection.cpp
 * \brief Object for holding a collection of ProtinMatches, PeptideMatches, and SpectrumMatches
 ********************************/
#include "ProteinMatchCollection.h"
#include "ProteinMatch.h"
#include "PeptideMatch.h"
#include "SpectrumMatch.h"

#include "Match.h"
#include "MatchCollection.h"

using namespace std;
using namespace Crux;

/**
 * \returns a blank ProteinMatchCollection
 */
ProteinMatchCollection::ProteinMatchCollection() {
}

/**
 * \returns a ProteinMatchCollection using a MatchCollection
 * TODO - remove this later
 */
ProteinMatchCollection::ProteinMatchCollection(
  MatchCollection* match_collection ///< matches to add
  ) {
  
  addMatches(match_collection);

}

  
/**
 * Default destructor
 */
ProteinMatchCollection::~ProteinMatchCollection() {

  // delete protein matches
  for (ProteinMatchIterator iter = proteinMatchBegin();
       iter != proteinMatchEnd();
       ++iter) {
    delete *iter;
  }

  // delete modified aa strings + peptide matches
  for (map<MODIFIED_AA_T*, PeptideMatch*, cmpSeq>::iterator iter = peptide_match_map_.begin();
       iter != peptide_match_map_.end();
       ++iter) {
      free(iter->first);
      delete iter->second;
  }

  // delete spectrum matches
  for (SpectrumMatchIterator iter = spectrumMatchBegin();
       iter != spectrumMatchEnd();
       ++iter) {
    delete *iter;
  }

}

/**
 * \returns the begin iterator for the ProteinMatch objects
 */
ProteinMatchIterator ProteinMatchCollection::proteinMatchBegin() {

  return protein_matches_.begin();
}

/**
 * \returns the end iterator for the ProteinMatch objects
 */
ProteinMatchIterator ProteinMatchCollection::proteinMatchEnd() {

  return protein_matches_.end();
}

/**
 * \returns the begin iterator for the PeptideMatch objects
 */
PeptideMatchIterator ProteinMatchCollection::peptideMatchBegin() {
  return peptide_matches_.begin();
}
  
/**
 * \returns the end iterator for the PeptideMatch objects
 */
PeptideMatchIterator ProteinMatchCollection::peptideMatchEnd() {
  return peptide_matches_.end();
}

/**
 * \returns the begin iterator for all of the SpectrumMatch objects
 */
SpectrumMatchIterator ProteinMatchCollection::spectrumMatchBegin() {
  return spectrum_matches_.begin();
}

/**
 * \returns the end iterator for all of the SpectrumMatch objects
 */
SpectrumMatchIterator ProteinMatchCollection::spectrumMatchEnd() {
  return spectrum_matches_.end();
}

/**
 * \returns the ProteinMatch for a Protein, creates a new
 * one if it is not found and create is true
 */  
ProteinMatch* ProteinMatchCollection::getProteinMatch(
  Protein* protein,  ///< Protein for which to find the protein match
  bool create ///< Create the ProteinMatch if it doesn't exist
  ) {

  string id = protein->getIdPointer();

  ProteinMatch* ans = getProteinMatch(id);
  if (ans == NULL) {
    if (create) {
      ans = new ProteinMatch(protein);
      protein_matches_.push_back(ans);
      protein_match_map_.insert(pair<string, ProteinMatch*>(id, ans));
    } else {
      carp(CARP_WARNING, "ProteinMatch for %s not found!", protein->getIdPointer());
    }
  }
  return ans;
}

/**
 * \returns the PeptideMatch for a peptide object, creating
 * if it doesn't exist and create is true
 */
PeptideMatch* ProteinMatchCollection::getPeptideMatch(
  Peptide* peptide, ///< peptide to find
  bool create ///< create if it doesn't exist?
  ) {


  MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();

  PeptideMatch* ans = getPeptideMatch(mod_seq);

  if (ans == NULL) {
    if (create) {
      ans = new PeptideMatch(peptide);
      peptide_matches_.push_back(ans);
      peptide_match_map_.insert(pair<MODIFIED_AA_T*, PeptideMatch*>(mod_seq, ans));
    } else {
      free(mod_seq);
      carp(CARP_FATAL, "Could not find peptidematch for sequence %s",
        peptide->getSequence());
    }
  } else {
    free(mod_seq);
  }

  return ans;
}

/**
 * \returns the ProteinMatch for a Protein, null if it doesn't exist
 */
ProteinMatch* ProteinMatchCollection::getProteinMatch(
  const string& id ///< id of the protein
  ) {
  map<string, ProteinMatch*>::iterator iter = protein_match_map_.find(id);
  return (iter != protein_match_map_.end()) ? iter->second : NULL;
}

/**
 * \returns the PeptideMatch for the sequence, null if it doesn't exist.
 */
PeptideMatch* ProteinMatchCollection::getPeptideMatch(
  MODIFIED_AA_T* mod_seq ///< Modified Sequence to find
  ) {
  map<MODIFIED_AA_T*, PeptideMatch*, cmpSeq>::iterator iter = peptide_match_map_.find(mod_seq);
  return (iter != peptide_match_map_.end()) ? iter->second : NULL;
}



/**
 * Helper method that adds a Crux match to the ProteinCollection,
 * creating the SpectrumMatch, PeptideMatch, and ProteinMatch objects.
 */
void ProteinMatchCollection::addMatch(
  MatchCollection* match_collection, ///< Collection from where the match came from
  Match* match ///< Match to add
  ){

  //create a spectrum match.
  Spectrum* spectrum = match->getSpectrum();
  SpectrumZState z_state = match->getZState();
  SpectrumMatch* spectrum_match = new SpectrumMatch(spectrum);
  spectrum_match->setFileIndex(match->getFileIndex());
  spectrum_match->setZState(z_state);
  spectrum_match->setScore(DELTA_CN, match->getDeltaCn());
  spectrum_match->setScore(BY_IONS_MATCHED, match->getBYIonMatched());
  spectrum_match->setScore(BY_IONS_TOTAL, match->getBYIonPossible());
  spectrum_matches_.push_back(spectrum_match);
  
  pair<int, int> scan_charge = make_pair(spectrum->getFirstScan(), z_state.getCharge());

  if (match->getLnExperimentSize() >= 0 &&
      spectrum_counts_.find(scan_charge) == spectrum_counts_.end()) {
    spectrum_counts_[scan_charge] = floor(exp(match->getLnExperimentSize()) + 0.5);
  }

  Peptide* peptide = match->getPeptide();
  PeptideMatch* peptide_match = getPeptideMatch(peptide);
  
  //add the spectrum match.
  peptide_match->addSpectrumMatch(spectrum_match);

  for (int score_idx = (int)SP;
       score_idx < (int)NUMBER_SCORER_TYPES;
       score_idx++) {

    SCORER_TYPE_T score_type = (SCORER_TYPE_T)score_idx;
    if (match_collection->getScoredType(score_type)) {
      peptide_match->setScore(score_type, match->getScore(score_type));
      peptide_match->setRank(score_type, match->getRank(score_type));
      spectrum_match->setScore(score_type, match->getScore(score_type));
      spectrum_match->setRank(score_type, match->getRank(score_type));
    }
  }

  for (PeptideSrcIterator src_iter = peptide->getPeptideSrcBegin();
       src_iter != peptide->getPeptideSrcEnd();
       ++src_iter) {
    PeptideSrc* src = *src_iter;
    Protein* protein = src->getParentProtein();
    ProteinMatch* protein_match = getProteinMatch(protein);
    protein_match->addPeptideMatch(peptide_match);
    peptide_match->addProteinMatch(protein_match, src);
  }
}

/**
 * Helper method that adds a Crux match to the ProteinMatchCollection,
 * Adding all SpectrumMatch, PeptideMatch, and ProteinMatches
 */
void ProteinMatchCollection::addMatches(
  MatchCollection* match_collection ///< collection to add
  ) {

  carp(CARP_DEBUG, "Adding %d matches to ProteinMatchCollection",
       match_collection->getMatchTotal());
  distinct_matches_ = match_collection->getHasDistinctMatches();

  MatchIterator match_iter(match_collection);
  while(match_iter.hasNext()) {
    addMatch(match_collection, match_iter.next());
  }

}

/**
 * Get the matches/spectrum as a map, where the key is <scan, chage>
 */
const map<pair<int, int>, int>& ProteinMatchCollection::getMatchesSpectrum() {
  return spectrum_counts_;
}


/**                                                                                                                                                                                                      
 * \returns whether matches are distinct are not                                                                                                                                                         
 */
bool ProteinMatchCollection::hasDistinctMatches() {
  return distinct_matches_;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
