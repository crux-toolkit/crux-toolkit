#include "PeptideMatch.h"
#include "SpectrumMatch.h"

using namespace Crux;

PeptideMatch::PeptideMatch() {

}

PeptideMatch::PeptideMatch(
  Peptide* peptide
  ) {

  setPeptide(peptide);

}

PeptideMatch::~PeptideMatch() {

}

void PeptideMatch::setPeptide(
  Peptide* peptide
  ) {

  peptide_ = peptide;

}

Peptide* PeptideMatch::getPeptide() {

  return(peptide_);
}

void PeptideMatch::addProteinMatch(
  ProteinMatch* protein_match,
  PeptideSrc* src
  ) {
  protein_matches_.push_back(protein_match);
  protein_match_to_peptide_src_[protein_match] = src;
}

PeptideSrc* PeptideMatch::getSrc(
  ProteinMatch* protein_match
  ) {

  return protein_match_to_peptide_src_[protein_match];

}

void PeptideMatch::addSpectrumMatch(
    SpectrumMatch* spectrum_match
    ) {
  spectrum_matches_.push_back(spectrum_match);
  spectrum_match->setPeptideMatch(this);  
}


SpectrumMatchIterator PeptideMatch::spectrumMatchBegin() {
  return spectrum_matches_.begin();
}
  
SpectrumMatchIterator PeptideMatch::spectrumMatchEnd() {
  return spectrum_matches_.end();
}

ProteinMatchIterator PeptideMatch::proteinMatchBegin() {
  return protein_matches_.begin();
}

ProteinMatchIterator PeptideMatch::proteinMatchEnd() {
  return protein_matches_.end();
}
#include "PeptideMatch.h"
#include "SpectrumMatch.h"

using namespace Crux;

PeptideMatch::PeptideMatch() {

}

PeptideMatch::PeptideMatch(
  Peptide* peptide
  ) {

  setPeptide(peptide);

}

PeptideMatch::~PeptideMatch() {

}

void PeptideMatch::setPeptide(
  Peptide* peptide
  ) {

  peptide_ = peptide;

}

Peptide* PeptideMatch::getPeptide() {

  return(peptide_);
}

void PeptideMatch::addProteinMatch(
  ProteinMatch* protein_match,
  PeptideSrc* src
  ) {
  protein_matches_.push_back(protein_match);
  protein_match_to_peptide_src_[protein_match] = src;
}

PeptideSrc* PeptideMatch::getSrc(
  ProteinMatch* protein_match
  ) {

  return protein_match_to_peptide_src_[protein_match];

}

void PeptideMatch::addSpectrumMatch(
    SpectrumMatch* spectrum_match
    ) {
  spectrum_matches_.push_back(spectrum_match);
  spectrum_match->setPeptideMatch(this);  
}


SpectrumMatchIterator PeptideMatch::spectrumMatchBegin() {
  return spectrum_matches_.begin();
}
  
SpectrumMatchIterator PeptideMatch::spectrumMatchEnd() {
  return spectrum_matches_.end();
}

ProteinMatchIterator PeptideMatch::proteinMatchBegin() {
  return protein_matches_.begin();
}

ProteinMatchIterator PeptideMatch::proteinMatchEnd() {
  return protein_matches_.end();
}
