/**
 * \file spectrum_collection.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * $Revision: 1.28 $
 * \brief Object for representing many spectra.
 *****************************************************************************/
#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include <stdio.h>
#include "objects.h"
#include "spectrum.h"
#include "carp.h"


/**
 * \returns An (empty) spectrum_collection object.
 */
SPECTRUM_COLLECTION_T* allocate_spectrum_collection(void);

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does *NOT* parse all spectra in the file. 
 * This will be done lazily depending on the subsequent method
 * calls (parse_spectrum_collection get_spectrum_collection_spectrum).
 * \returns  SPECTRUM_COLLECTION_T
 */
SPECTRUM_COLLECTION_T* new_spectrum_collection(
  const char* filename///< The spectrum collection filename. -in
  );

/**
 * Frees an allocated spectrum_collection object.
 */
void free_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum collection to free - in
);
/**
 * Prints a spectrum_collection object to file.
 */
void print_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< spectrum_collection to print -in 
  FILE* file ///< file for output -out
  );

/**
 * Copies spectrum_collection object from src to dest.
 *  must pass in a memory allocated SPECTRUM_COLLECTION_T* dest
 */
void copy_spectrum_collection(
  SPECTRUM_COLLECTION_T* src,///< spectrum to copy from -in
  SPECTRUM_COLLECTION_T* dest///< spectrum to copy to -out
  );

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
BOOLEAN_T parse_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< empty spectrum to parse into -out
);

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
BOOLEAN_T get_spectrum_collection_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< The spectrum collection -out
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  SPECTRUM_T* spectrum ///< The (empty) allocated SPECTRUM_T object -in
  );

/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum in correct order into the spectra array
 * spectrum must be heap allocated
 *\returns TRUE if succeed to add, else FALSE 
 */
BOOLEAN_T add_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to add to spectrum_collection -in
); 

/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum to the end of the spectra array
 * should only be used when the adding in increasing scan num order
 * when adding in random order should use add_spectrum
 * spectrum must be heap allocated
 *\returns TRUE if succeed to add, else FALSE 
 */
BOOLEAN_T add_spectrum_to_end(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to add to spectrum_collection -in
  );

/**
 * Removes a spectrum from the spectrum_collection.
 */
void remove_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to be removed from spectrum_collection -in
  ); 

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**  ////// TESTME////
 * \sets the filename of the ms2 file the spectra were parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_collection_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_collection_new_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save filename -out
  char* filename ///< filename -in
  );

/**
 * \returns the filename of the ms2 file the spectra were parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
void set_spectrum_collection_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save filename -out
  char* filename ///< filename -in
  );

/**
 * \returns the filename of the ms2 file the spectra was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* get_spectrum_collection_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum collection's filename -in 
  );

/**
 * \returns the current number of spectrum in the spectrum_collection
 */
int get_spectrum_collection_num_spectra(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
  );

/**
 * \returns The current number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra) in the spectrum_collection
 */
int get_spectrum_collection_num_charged_spectra(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection save filename -in
  );

/**
 * \returns the comments from the spectrum_collection
 * the return char* points to a newly heap allocated copy of the comments
 * user must free the new string object
 */
char* get_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
  );


/**
 * sets the comment of the spectrum_collection, If comments exist add this to the end
 * copies the new_comment into a newly heap allocated copy of the comment
 */
void set_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save comment -in                                         
  char* new_comment ///< the new comments to be copied
  );

/**
 * \returns TRUE if the spectrum_collection file has been parsed
 */
BOOLEAN_T get_spectrum_collection_is_parsed(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in
);

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
FILE** get_spectrum_collection_psm_result_filenames(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection -in
  char* psm_result_folder_name, ///< the folder name for where the result file should be placed -in
  char*** psm_result_filenames, ///< pointer to be set to the array of filenames for both the target and decoy psm results -out
  int number_decoy_set,  ///< the number of decoy sets to produce -in
  char* file_extension ///< the file extension of the spectrum file(i.e. ".ms2") -in
  );


/******************************************************************************/

/**
 * Instantiates a new spectrum_iterator object from spectrum_collection.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
SPECTRUM_ITERATOR_T* new_spectrum_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectrum_collection to iterate -in
);        

/**
 * Instantiates a new spectrum_iterator object from
 * spectrum_collection.  This iterator returns unique spectrum-charge
 * pairs (e.g.a spectrum to be searched as +2 and +3 is returned once
 * as +2 and once as +3).  The charge is returned by setting the int
 * pointer in the argument list.  The iterator also filters spectra by
 * mass so that none outside the spectrum-min-mass--spectrum-max-mass
 * range (as defined in parameter.c).
 * \returns a SPECTRUM_ITERATOR_T object.
 */
FILTERED_SPECTRUM_CHARGE_ITERATOR_T* new_filtered_spectrum_charge_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectra to iterate over
);        

/**
 * Frees an allocated spectrum_iterator object.
 */
void free_spectrum_iterator(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< free spectrum_iterator -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T spectrum_iterator_has_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< is there a next spectrum? -in
);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T filtered_spectrum_charge_iterator_has_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator);

/**
 * The basic iterator function next.
 */
SPECTRUM_T* spectrum_iterator_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< return the next spectrum -in
);

/**
 * The basic iterator function next.  Also returns the charge state to
 * use for this spectrum.
 */
SPECTRUM_T* filtered_spectrum_charge_iterator_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator,///< return spec from here -in
  int* charge);                 ///< put charge here -out


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
