/**
 * \file ProteinMatch.cpp
 * \brief Object for holding protein scores
 ********************************/
#include "ProteinMatch.h"

#include "Protein.h"

using namespace Crux;
using namespace std;

/**
 * \returns an empty ProteinMatch object
 */
ProteinMatch::ProteinMatch() {
  ;
}
  
/**
 * \returns a ProteinMatch object with the protein set
 */
ProteinMatch::ProteinMatch(
  Protein* protein ///< Protein to set
  ) {
  protein_ = protein;
}
    
/**
 * Destructor
 */
ProteinMatch::~ProteinMatch() {

}

/** 
 * \returns the protein for the match
 */
Crux::Protein* ProteinMatch::getProtein() {
  return protein_;
}

/**
 * sets the protein for the match
 */
void ProteinMatch::setProtein(
  Crux::Protein* protein ///< Protein to set
) {
  protein_ = protein;
}

/**
 * \returns the begining of the peptidematches vector
 */
PeptideMatchIterator ProteinMatch::peptideMatchBegin() {
  return peptide_matches_.begin();
}

/**
 * \returns the end of the peptide matches vector
 */
PeptideMatchIterator ProteinMatch::peptideMatchEnd() {
  return peptide_matches_.end();
}

/**
 * adds a peptide match to the protein match
 */
void ProteinMatch::addPeptideMatch(
  PeptideMatch* peptide_match ///<PeptideMatch to add
  ) {
  for (deque<PeptideMatch*>::iterator iter = peptide_matches_.begin();
       iter != peptide_matches_.end();
       ++iter) {
    if (*iter == peptide_match) {
      return;
    }
  }
  peptide_matches_.push_back(peptide_match);
}

/**
 * \returns the id for this protein
 */
string ProteinMatch::getId() const {
  return string(protein_->getIdPointer());
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
