/**
 * \file PWIZSpectrumCollection.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 21 April 2011
 * \brief Class for accessing spectra using the proteoewizard library.
 */
#ifndef PWIZ_SPECTRUM_COLLECTION_H
#define PWIZ_SPECTRUM_COLLECTION_H

#include "SpectrumCollection.h"

#include "pwiz/data/msdata/MSDataFile.hpp"

/**
 * \class SpectrumCollection
 * \brief An abstract class for accessing spectra from a file.
 */
class PWIZSpectrumCollection : public Crux::SpectrumCollection {

 protected:
  pwiz::msdata::MSDataFile* reader_;
  
  /**
   * Parses the first/last scan from the title
   * \returns whether parse was successful.
   * For MGF files that place their scan numbers in the title string
   */
  bool parseFirstLastScanFromTitle(
    std::string& scan_title_str, ///< title string to parse -in
    int& first_scan, ///< first scan -out
    int& last_scan ///< last scan -out
  );
  
 public:
  /**
   * Constructor sets filename and initializes member variables.
   */
  PWIZSpectrumCollection(
    const char* filename ///< The spectrum collection filename. -in
  );

  /**
   * Default Constructor
   */
  virtual ~PWIZSpectrumCollection();

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
