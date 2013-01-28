#include "SpectrumMatch.h"
#include "SpectrumZState.h"

using namespace Crux;

SpectrumMatch::SpectrumMatch() {
}

SpectrumMatch::SpectrumMatch(
  Spectrum* spectrum
  ) {

  setSpectrum(spectrum);

}
  
SpectrumMatch::~SpectrumMatch() {
}
  
void SpectrumMatch::setPeptideMatch(
  PeptideMatch* peptide_match
  ) {
  peptide_match_ = peptide_match;
}

PeptideMatch* SpectrumMatch::getPeptideMatch() {
  return peptide_match_;
}


void SpectrumMatch::setSpectrum(
  Spectrum* spectrum
  ) {

  spectrum_ = spectrum;
}

Spectrum* SpectrumMatch::getSpectrum() {
  return spectrum_;
}

void SpectrumMatch::setZState(
  SpectrumZState& zstate
  ) {

  zstate_ = zstate;
}

SpectrumZState& SpectrumMatch::getZState() {

  return zstate_;
}
#include "SpectrumMatch.h"
#include "SpectrumZState.h"

using namespace Crux;

SpectrumMatch::SpectrumMatch() {
}

SpectrumMatch::SpectrumMatch(
  Spectrum* spectrum
  ) {

  setSpectrum(spectrum);

}
  
SpectrumMatch::~SpectrumMatch() {
}
  
void SpectrumMatch::setPeptideMatch(
  PeptideMatch* peptide_match
  ) {
  peptide_match_ = peptide_match;
}

PeptideMatch* SpectrumMatch::getPeptideMatch() {
  return peptide_match_;
}


void SpectrumMatch::setSpectrum(
  Spectrum* spectrum
  ) {

  spectrum_ = spectrum;
}

Spectrum* SpectrumMatch::getSpectrum() {
  return spectrum_;
}

void SpectrumMatch::setZState(
  SpectrumZState& zstate
  ) {

  zstate_ = zstate;
}

SpectrumZState& SpectrumMatch::getZState() {

  return zstate_;
}
