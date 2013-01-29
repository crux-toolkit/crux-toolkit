#include "ProteinMatch.h"

#include "Protein.h"

using namespace Crux;
using namespace std;

ProteinMatch::ProteinMatch() {
  ;
}

ProteinMatch::ProteinMatch(
  Protein* protein
) {
  protein_ = protein;
}
  
ProteinMatch::~ProteinMatch() {

}

PeptideMatchIterator ProteinMatch::peptideMatchBegin() {
  return peptide_matches_.begin();
}

PeptideMatchIterator ProteinMatch::peptideMatchEnd() {
  return peptide_matches_.end();
}

void ProteinMatch::addPeptideMatch(
  PeptideMatch* peptide_match
  ) {

  peptide_matches_.push_back(peptide_match);

}

string ProteinMatch::getId() {
  string ans = protein_->getIdPointer();
  return ans;

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
