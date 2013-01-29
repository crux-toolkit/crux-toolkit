/**
 * \file SpectrumMatch.h
 * $Revision: 1.00 $
 * \brief Object for holding spectrum scores
 ******************************************************/
#ifndef SPECTRUMMATCH_H_
#define SPECTRUMMATCH_H_

#include "match_objects.h"
#include "SpectrumZState.h"
#include "AbstractMatch.h"

class SpectrumMatch : public AbstractMatch {

 protected:
  PeptideMatch* peptide_match_;
  Crux::Spectrum* spectrum_;
  SpectrumZState zstate_;

 public:
  SpectrumMatch();
  SpectrumMatch(
    Crux::Spectrum* spectrum_
  );

  ~SpectrumMatch();

  void setPeptideMatch(
    PeptideMatch* peptide_match
  );

  PeptideMatch* getPeptideMatch();

  void setSpectrum(
    Crux::Spectrum* spectrum
  );

  Crux::Spectrum* getSpectrum();

  void setZState(
    SpectrumZState& zstate
  );

  SpectrumZState& getZState();

};

#endif //SPECTRUMMATCH_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
