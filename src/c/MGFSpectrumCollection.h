/**
 * \file MGFSpectrumCollection.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 April 2011
 * \brief Class for accessing spectra in .mgf (Mascot general format) files.
 */
#ifndef MGF_SPECTRUM_COLLECTION_H
#define MGF_SPECTRUM_COLLECTION_H

#include <stdio.h>
#include "objects.h"
#include "SpectrumCollection.h"

#include <vector>

/**
 * \class SpectrumCollection
 * \brief An abstract class for accessing spectra from a file.
 */
class MGFSpectrumCollection : public Crux::SpectrumCollection {

 protected:
  int cur_spec_number_;  ///< number spec sequentially in lieu of scan numbers

 public:
  /**
   * Constructor sets filename and initializes member variables.
   */
  MGFSpectrumCollection(
    const char* filename ///< The spectrum collection filename. -in
  );

  /**
   * Parses all the spectra from file designated by the filename member
   * variable.
   * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
   */
  virtual bool parse();

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.
   * \returns The newly allocated Spectrum or NULL if scan number not found.
   */
  virtual Crux::Spectrum* getSpectrum(
    int first_scan      ///< The first scan of the spectrum to retrieve -in
  );

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.  Removes any existing information in the
   * given spectrum.
   * \returns True if the spectrum was allocated, false on error.
   */
  virtual bool getSpectrum(
    int first_scan,      ///< The first scan of the spectrum to retrieve -in
    Crux::Spectrum* spectrum   ///< Put the spectrum info here
  );

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
