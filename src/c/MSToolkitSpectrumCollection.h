/**
 * \file MSToolkitSpectrumCollection.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 20 April 2011
 * \brief Class for accessing spectra using the MSToolkit library.
 */
#ifndef MSTOOLKIT_SPECTRUM_COLLECTION_H
#define MSTOOLKIT_SPECTRUM_COLLECTION_H

//#include <stdio.h>
//#include "objects.h"
#include "SpectrumCollection.h"

/**
 * \class SpectrumCollection
 * \brief An abstract class for accessing spectra from a file.
 */
class MSToolkitSpectrumCollection : public Crux::SpectrumCollection {

 protected:

 public:
  /**
   * Constructor sets filename and initializes member variables.
   */
  MSToolkitSpectrumCollection(
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
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
