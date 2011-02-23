/**
 * \file SpectrumCollection.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.28 $
 * \brief Object for representing many spectra.
 *****************************************************************************/
#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include <stdio.h>
#include "objects.h"
#include "Spectrum.h"
#include "carp.h"

#include <vector>


typedef std::vector<Spectrum*>::iterator SpectrumIterator; 

/**
 * \class spectrum_collection 
 * \brief A object to group together one or more spectrum objects.
 */
class SpectrumCollection {

 friend class FilteredSpectrumChargeIterator;

 protected:
  static const unsigned int MAX_COMMENT = 1000; ///< max length of comment
  std::vector<Spectrum*> spectra_;  ///< The spectrum peaks
  char* filename_;     ///< Optional filename
  char comment_[MAX_COMMENT];    ///< The spectrum_collection header lines
  int num_charged_spectra_;
  BOOLEAN_T is_parsed_; ///< Have we parsed all the spectra from the file? 
  
  /* Private functions */

  // TESTME might have to check for bad format
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
   * Default constructor.
   */
  SpectrumCollection();

  SpectrumCollection(
    const char* filename///< The spectrum collection filename. -in
  );

  /**
   * Copy constructor.
   */
  SpectrumCollection(SpectrumCollection& old_collection);

  /**
   * Default destructor.
   */
  ~SpectrumCollection();


  SpectrumIterator begin();
  SpectrumIterator end();
  

  /**
   * Prints a spectrum_collection object to file.
   */
  void printSpectrumCollection(
    FILE* file ///< file for output -out
    );

  /**
   * Parses all the spectra from file designated by the filename member
   * variable.
   * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
   */
  bool parse();

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.
   * \returns The newly allocated Spectrum or NULL if scan number not found.
   */
  Spectrum* getSpectrum(
    int first_scan      ///< The first scan of the spectrum to retrieve -in
  );

  /**
   * Parses a single spectrum from a spectrum_collection with first scan
   * number equal to first_scan.  Removes any existing information in the
   * given spectrum.
   * \returns True if the spectrum was allocated, false on error.
   */
  bool getSpectrum(
    int first_scan,      ///< The first scan of the spectrum to retrieve -in
    Spectrum* spectrum   ///< Put the spectrum info here
  );

  /**
   * Adds a spectrum to the spectrum_collection.
   * adds the spectrum in correct order into the spectra array
   * spectrum must be heap allocated
   *\returns TRUE if succeed to add, else FALSE 
   */
  bool addSpectrum(
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
  bool addSpectrumToEnd(
    Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  );

  /**
   * Removes a spectrum from the spectrum_collection.
   */
  void removeSpectrum(
    Spectrum* spectrum ///< spectrum to be removed from spectrum_collection -in
  ); 

  /**
   * Additional get and set methods
   */

  /**  ////// TESTME////
   * \sets the filename of the ms2 file the spectra were parsed
   * this function should be used only the first time the filename is set
   * to change existing filename use set_spectrum_collection_filename
   * copies the value from arguement char* filename into a heap allocated memory
   */
  void setNewFilename(
    char* filename ///< filename -in
  );

  /**
   * \returns the filename of the ms2 file the spectra were parsed
   * returns a char* to a heap allocated copy of the filename
   * user must free the memory
   */
  char* getFilename();

  void setFilename(
  char* filename ///< filename -in
  );

  /**
   * \returns the filename of the ms2 file the spectra was parsed
   * returns a char* to a heap allocated copy of the filename
   * user must free the memory
   */
  char* getSpectrumCollectionFilename();

  /**
   * \returns the current number of spectrum in the spectrum_collection
   */
  int getNumSpectra();

  /**
   * \returns The current number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra) in the spectrum_collection
   */
  int getNumChargedSpectra();

  /**
   * \returns the comments from the spectrum_collection
   * the return char* points to a newly heap allocated copy of the comments
   * user must free the new string object
   */
  char* getComment();

  /**
   * sets the comment of the spectrum_collection, If comments exist add this to the end
   * copies the new_comment into a newly heap allocated copy of the comment
   */
  void setComment(
    char* new_comment ///< the new comments to be copied
  );

  /**
   * \returns TRUE if the spectrum_collection file has been parsed
   */
  bool getIsParsed();

  /**
   * Takes the spectrum file name and creates a file with unique filenames.
   * The method will create one file for PSM result serializations for the target sequence
   * and #number_decoy_set number of files for decoy PSM result serialization.  
   * Thus, the FILE* array will contain,
   * at index 0, the target file and the fallowing indicies the decoy files.
   *
   * Template: "fileName_XXXXXX", where XXXXXX is random generated to be unique.
   * Also,sets psm_result_filenames pointer to the array of filenames for both the target and decoy psm results
   * The name array is heap allocated, thus user must free it. Size is number_decoy_set +1 (for target)
   *\returns file handler array to the newly created files(target& decoy) and sets psm_result_filename.
   */
  FILE** getPsmResultFilenames(
    char* psm_result_folder_name, ///< the folder name for where the result file should be placed -in
    char*** psm_result_filenames, ///< pointer to be set to the array of filenames for both the target and decoy psm results -out
    int number_decoy_set,  ///< the number of decoy sets to produce -in
    char* file_extension ///< the file extension of the spectrum file(i.e. ".ms2") -in
  );

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
