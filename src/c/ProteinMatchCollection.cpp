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
  Crux::Protein* protein,  ///< Protein for which to find the protein match
  bool create ///< Create the ProteinMatch if it doesn't exist
  ) {

  string id = protein->getIdPointer();

  ProteinMatch* ans = getProteinMatch(id);
  if (ans == NULL) {
    if (create) {
      ans = new ProteinMatch(protein);
      protein_matches_.push_back(ans);
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
    } else {
      carp(CARP_FATAL, 
        "Could not find peptidematch for sequence %s", 
	   peptide->getSequence());
    }
  }

  free(mod_seq);

  return ans;
}

/**
 * \returns the ProteinMatch for a Protein, null if it doesn't exist
 */
ProteinMatch* ProteinMatchCollection::getProteinMatch(
  const std::string& id ///< id of the protein
  ) {

  for (ProteinMatchIterator iter = proteinMatchBegin(); iter != proteinMatchEnd(); ++iter) {
    if ((*iter)->getId() == id) {
      return *iter;
    }
  }
  return NULL;
}

/**
 * \returns the PeptideMatch for the sequence, null if it doesn't exist.
 */
PeptideMatch* ProteinMatchCollection::getPeptideMatch(
  MODIFIED_AA_T* mod_seq ///< Modified Sequence to find
  ) {
  int ans = -1;
  for (size_t idx = 0; idx < peptide_matches_.size(); idx++) {
    MODIFIED_AA_T* cur_seq = peptide_matches_[idx]->getPeptide()->getModifiedAASequence();
    if (equal_seq(cur_seq, mod_seq)) {
      ans = idx;
      break;
    }
  }

  if (ans == -1) { 
    return NULL;
  } else {
    return peptide_matches_[ans];
  }


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
  SpectrumMatch* spectrum_match = new SpectrumMatch(match->getSpectrum());
  spectrum_match->setZState(match->getZState());
  spectrum_match->setScore(DELTA_CN, match->getDeltaCn());
  spectrum_match->setScore(BY_IONS_MATCHED, match->getBYIonMatched());
  spectrum_match->setScore(BY_IONS_TOTAL, match->getBYIonPossible());
  spectrum_matches_.push_back(spectrum_match);

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

  MatchIterator match_iter(match_collection);
  while(match_iter.hasNext()) {
    addMatch(match_collection, match_iter.next());
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
