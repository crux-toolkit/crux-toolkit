#include "ProteinMatchCollection.h"
#include "ProteinMatch.h"
#include "PeptideMatch.h"
#include "SpectrumMatch.h"

#include "Match.h"
#include "MatchCollection.h"

using namespace Crux;

ProteinMatchCollection::ProteinMatchCollection() {
}

ProteinMatchCollection::ProteinMatchCollection(
  MatchCollection* match_collection
  ) {

  addMatches(match_collection);

}

ProteinMatchCollection::~ProteinMatchCollection() {
}

ProteinMatchIterator ProteinMatchCollection::proteinMatchBegin() {

  return protein_matches_.begin();
}

ProteinMatchIterator ProteinMatchCollection::proteinMatchEnd() {

  return protein_matches_.end();
}

PeptideMatchIterator ProteinMatchCollection::peptideMatchBegin() {
  return peptide_matches_.begin();
}

PeptideMatchIterator ProteinMatchCollection::peptideMatchEnd() {
  return peptide_matches_.end();
}

SpectrumMatchIterator ProteinMatchCollection::spectrumMatchBegin() {
  return spectrum_matches_.begin();
}

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

  char* seq = peptide->getSequence();
  string seq_str = seq;
  free(seq);

  PeptideMatch* ans = getPeptideMatch(seq_str);

  if (ans == NULL) {
    if (create) {
      ans = new PeptideMatch(peptide);
      peptide_matches_.push_back(ans);
    } else {
      carp(CARP_WARNING, 
        "Could not find peptidematch for sequence %s", 
        seq_str.c_str());
    }
  }

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
  const std::string& sequence ///<sequence to find
  ) {

  for (size_t idx = 0; idx < peptide_matches_.size();idx++) {
    char* seq = peptide_matches_[idx]->getPeptide()->getSequence();
    string seq_str = seq;
    free(seq);

    if (seq_str == sequence) {
      return peptide_matches_[idx];
    }
  }
  return NULL;
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
      spectrum_match->setScore(score_type, match->getScore(score_type));
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
