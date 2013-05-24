/**
 * \file PeptideMatch.cpp
 * \brief Object for holding peptide scores
 ********************************/
#include "PeptideMatch.h"
#include "SpectrumMatch.h"

using namespace Crux;
using namespace std;

/**
 * \returns an empty PeptideMatch
 */
PeptideMatch::PeptideMatch() {

}

/**
 * \returns a PeptideMatch with the peptide set
 */
PeptideMatch::PeptideMatch(
  Peptide* peptide  ///< Peptide to set
  ) {

  setPeptide(peptide);

}

/**
 * Default destructor
 */
PeptideMatch::~PeptideMatch() {

}

/**
 * sets the peptide for this match
 */
void PeptideMatch::setPeptide(
  Peptide* peptide ///< peptide to set
  ) {

  peptide_ = peptide;

}

/**
 * \returns the peptide for the match
 */
Peptide* PeptideMatch::getPeptide() {

  return(peptide_);
}

/**
 * adds a ProteinMatch to this PeptideMatch
 */
void PeptideMatch::addProteinMatch(
  ProteinMatch* protein_match, ///< ProteinMatch to set
  PeptideSrc* src ///< PeptideSrc 
  ) {
  for (deque<ProteinMatch*>::iterator iter = protein_matches_.begin();
       iter != protein_matches_.end();
       ++iter) {
    if (*iter == protein_match) {
      return;
    }
  }
  protein_matches_.push_back(protein_match);
  protein_match_to_peptide_src_[protein_match] = src;
}

/**
 * \returns the associated PeptideSrc for the ProteinMatch
 */
PeptideSrc* PeptideMatch::getSrc(
  ProteinMatch* protein_match ///< gets the associated PeptideSrc
  ) {

  return protein_match_to_peptide_src_[protein_match];

}

/**
 * Sets the spectrum match
 */
void PeptideMatch::addSpectrumMatch(
  SpectrumMatch* spectrum_match ///< SpectrumMatch to add
  ) {
  for (deque<SpectrumMatch*>::iterator iter = spectrum_matches_.begin();
       iter != spectrum_matches_.end();
       ++iter) {
    if (*iter == spectrum_match) {
      return;
    }
  }
  spectrum_matches_.push_back(spectrum_match);
  spectrum_match->setPeptideMatch(this);  
}

  
/**
 * \returns the begin iterator of the spectrummatches
 */
SpectrumMatchIterator PeptideMatch::spectrumMatchBegin() {
  return spectrum_matches_.begin();
}

/**
 * \returns the begin iterator of the spectrum matches
 */
SpectrumMatchIterator PeptideMatch::spectrumMatchEnd() {
  return spectrum_matches_.end();
}

/**
 * \returns the begin iterator of the spectrum matches
 */
ProteinMatchIterator PeptideMatch::proteinMatchBegin() {
  return protein_matches_.begin();
}

/**
 * \returns the end iterator of the spectrum matches
 */
ProteinMatchIterator PeptideMatch::proteinMatchEnd() {
  return protein_matches_.end();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
