/**
 * \file SpectrumCollection.h 
 * AUTHOR: Barbara Frewen
 * CREATE DATE: 14 April 2011
 * \brief An abstract class for accessing spectra from a file.
 */
#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include <stdio.h>
#include "objects.h"
#include "Spectrum.h"

#include <vector>


typedef std::vector<Spectrum*>::iterator SpectrumIterator; 

/**
 * \class SpectrumCollection
 * \brief An abstract class for accessing spectra from a file.
 */
class SpectrumCollection {

 friend class FilteredSpectrumChargeIterator;

 protected:
  std::vector<Spectrum*> spectra_;  ///< spectra from the file
  std::string filename_;                  ///< filename
  bool is_parsed_;      ///< file has been read and spectra_ populated 
  int num_charged_spectra_;  ///< sum of all charge states from all spectra
  
  /**
   * Base class constructor is protected.  Sets filename and
   * initializes member variables.
   */
  SpectrumCollection(
    const char* filename///< The spectrum collection filename. -in
  );

   /**
   * Adds a spectrum to the spectrum_collection.
   * adds the spectrum in correct order into the spectra array
   * spectrum must be heap allocated
   *\returns TRUE if succeed to add, else FALSE 
   */
  void addSpectrum(
    Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  ); 

  /**
   * Adds a spectrum to the spectrum_collection.
   * adds the spectrum to the end of the spectra array
   * should only be used when the adding in increasing scan num order
   * when adding in random order should use add_spectrum
   * spectrum must be heap allocated
   *\returns TRUE if succeed to add, else FALSE 
   */
  void addSpectrumToEnd(
    Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  );

  /**
   * Removes a spectrum from the spectrum_collection.
   */
  void removeSpectrum(
    Spectrum* spectrum ///< spectrum to be removed from spectrum_collection -in
  ); 

 public:

  /**
   * Copy constructor.  Deep copy.
   */
  SpectrumCollection(SpectrumCollection& old_collection);
  /**
   * Default destructor.
   */
  virtual ~SpectrumCollection();

 
  SpectrumIterator begin();
  SpectrumIterator end();
  

  /**
   * Parses all the spectra from file designated by the filename member
   * variable.
   * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
   */
  virtual bool parse() = 0;

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.
   * \returns The newly allocated Spectrum or NULL if scan number not found.
   */
  virtual Spectrum* getSpectrum(
    int first_scan      ///< The first scan of the spectrum to retrieve -in
  ) = 0;

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.  Removes any existing information in the
   * given spectrum.
   * \returns True if the spectrum was allocated, false on error.
   */
  virtual bool getSpectrum(
    int first_scan,      ///< The first scan of the spectrum to retrieve -in
    Spectrum* spectrum   ///< Put the spectrum info here
  ) = 0;

  /**
   * \returns A pointer to the name of the file containing these spectra.
   */
  const char* getFilename();

  /**
   * \returns The current number of spectrum in the
   * spectrum_collection.  Assumes it has been parsed.
   */
  int getNumSpectra();

  /**
   * \returns The current number of spectra assuming differnt
   * charge(i.e. one spectrum with two charge states are counted as
   * two spectra) in the spectrum_collection.  Assumes file has been
   * parsed.  
   */
  int getNumChargedSpectra();

  /**
   * \returns TRUE if the spectrum_collection file has been parsed
   */
  bool getIsParsed();
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
