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
  PeptideMatch* peptide_match_; ///< PeptideMatch which this SpectrumMatch matches to
  Crux::Spectrum* spectrum_; ///< Spectrum object
  SpectrumZState zstate_; ///< Current zstate
  int file_idx_; /// file index
 public:

  /**
   * \returns an empty SpectrumMatch
   */
  SpectrumMatch();

  /**
   * \returns a SpectrumMatch, setting the spectrum
   */
  SpectrumMatch(
    Crux::Spectrum* spectrum_ ///< Spectrum object for this match
  );

  /**
   * Destructor
   */
  virtual ~SpectrumMatch();

  /**
   * sets the peptide match for this spectrummatch
   */
  void setPeptideMatch(
    PeptideMatch* peptide_match ///< peptide match to set
  );

  /**
   * \returns the associated PeptideMatch
   */
  PeptideMatch* getPeptideMatch();

  /**
   * sets the spectrum object
   */
  void setSpectrum(
    Crux::Spectrum* spectrum ///< spectrum to set
  );

  /**
   * \returns the spectrum
   */
  Crux::Spectrum* getSpectrum();

  /**
   * sets the zstate
   */
  void setZState(
    SpectrumZState& zstate ///< zstate to set
  );

  /**
   * \returns the zstate
   */
  SpectrumZState& getZState();

  /**
   * sets the file index for this spectrum match
   */
  void setFileIndex(int file_idx);
  
  /**
   * \returns the file index for this spectrum match
   */
  int getFileIndex();

  /**
   *\returns the file path for this spectrum match
   */
  std::string getFilePath();
  
  
  
};

#endif //SPECTRUMMATCH_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
