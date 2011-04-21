/**
 * \file MS2SpectrumCollection.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.28 $
 * \brief Object for Reading spectra from a .ms2 file.
 */
#ifndef MS2_SPECTRUM_COLLECTION_H
#define MS2_SPECTRUM_COLLECTION_H

#include "SpectrumCollection.h"


typedef std::vector<Spectrum*>::iterator SpectrumIterator; 

/**
 * \class spectrum_collection 
 * \brief A object to group together one or more spectrum objects.
 */
class MS2SpectrumCollection : public SpectrumCollection {


 protected:
  static const unsigned int MAX_COMMENT = 1000; ///< max length of comment
  char comment_[MAX_COMMENT];    ///< The spectrum_collection header lines
  
  /**
   * parses 'H' line into the spectrum_collection comments
   * all reminding comments are ignored if max length of comments are reached
   */
  void parseHeaderLine(
    FILE* file
  );

  static long binarySearchSpectrum(
    FILE* file, 
    int first_scan
  );

  static int matchFirstScanLine(
    char* line,
    int buf_length,
    int query_first_scan
  );

 public:

  /**
   * Instantiates a new spectrum_collection object from a filename. 
   * Does not parse all spectra in the file. 
   * This will be done lazily depending on the subsequent method
   * calls (parse() or  getSpectrum()).
   */
  MS2SpectrumCollection(
    const char* filename ///< The spectrum collection filename.
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
  virtual Spectrum* getSpectrum(
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
    Spectrum* spectrum   ///< Put the spectrum info here
  );
  

  /**
   * \returns The comments (header lines) from the ms2 file.
   */
  const char* getComment();

  /**
   * Sets or adds to the the comment of the ms2 file.
   */
  void setComment(
    const char* new_comment ///< the new comments to be copied
  );


}; // class

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
