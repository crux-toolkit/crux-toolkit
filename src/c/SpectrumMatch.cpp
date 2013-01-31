/**
 * \file SpectrumMatch.cpp
 * \brief Object for holding spectrum scores
 *******************************************/
#include "SpectrumMatch.h"
#include "SpectrumZState.h"

using namespace Crux;

/**
 * \returns an empty SpectrumMatch
 */
SpectrumMatch::SpectrumMatch() {
}

/**
 * \returns a SpectrumMatch, setting the spectrum
 */
SpectrumMatch::SpectrumMatch(
  Spectrum* spectrum /// < Spectrum object for this match
  ) {

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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
