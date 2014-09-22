/**
 * \file SpectrumMatch.cpp
 * \brief Object for holding spectrum scores
 *******************************************/
#include "SpectrumMatch.h"
#include "SpectrumZState.h"

#include "Match.h"
using namespace Crux;

/**
 * \returns an empty SpectrumMatch
 */
SpectrumMatch::SpectrumMatch():
  file_idx_(-1) {
}

/**
 * \returns a SpectrumMatch, setting the spectrum
 */
SpectrumMatch::SpectrumMatch(
  Spectrum* spectrum /// < Spectrum object for this match
): file_idx_(-1) {

  setSpectrum(spectrum);

}

/**
 * Destructor
 */  
SpectrumMatch::~SpectrumMatch() {
}
  
/**
 * sets the peptide match for this spectrummatch
 */
void SpectrumMatch::setPeptideMatch(
  PeptideMatch* peptide_match ///< peptide match to set
  ) {
  peptide_match_ = peptide_match;
}

/**
 * \returns the associated PeptideMatch
 */
PeptideMatch* SpectrumMatch::getPeptideMatch() {
  return peptide_match_;
}

/**
 * sets the spectrum object
 */
void SpectrumMatch::setSpectrum(
  Spectrum* spectrum ///< spectrum to set
  ) {

  spectrum_ = spectrum;
}

/** 
 * \returns the Spectrum
 */
Spectrum* SpectrumMatch::getSpectrum() {
  return spectrum_;
}

/**
 * sets the zstate
 */
void SpectrumMatch::setZState(
  SpectrumZState& zstate ///< zstate to set
  ) {

  zstate_ = zstate;
}

/**
 * \returns the zstate
 */
SpectrumZState& SpectrumMatch::getZState() {

  return zstate_;
}

void SpectrumMatch::setFileIndex(int file_idx) {
  file_idx_ = file_idx;
}

int SpectrumMatch::getFileIndex() {
  return(file_idx_);
}

std::string SpectrumMatch::getFilePath() {
  return(Match::getFilePath(file_idx_));
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
