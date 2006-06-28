/**
 * \file spectrum_collection.h 
 * $Revision: 1.7 $
 * \brief Object for representing many spectra.
 *****************************************************************************/
#ifndef SPECTRUM_COLLECTION_H
#define SPECTRUM_COLLECTION_H

#include "spectrum.h"
#include <stdio.h>

/**
 * \typedef SPECTRUM_COLLECTION_T 
 */
typedef struct spectrum_collection SPECTRUM_COLLECTION_T;

/**
 * \typedef SPECTRUM_ITERATOR_T An object to iterate over the spectra in a
 * spectrum_collection object.
 */
typedef struct spectrum_iterator SPECTRUM_ITERATOR_T;

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
  SPECTRUM_COLLECTION_T* 
    spectrum_collection,  ///< An (empty) allocated spectrum_collection
  char* filename);        ///< The spectrum collection filename.

/**
 * Frees an allocated spectrum_collection object.
 */
void free_spectrum_collection(SPECTRUM_COLLECTION_T* spectrum_collection);

/**
 * Prints a spectrum_collection object to file.
 */
void print_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection, 
  FILE* file);

/**
 * Copies spectrum_collection object from src to dest.
 */
void copy_spectrum_collection(
  SPECTRUM_COLLECTION_T* src,
  SPECTRUM_COLLECTION_T* dest);

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
BOOLEAN_T parse_spectrum_collection(
    SPECTRUM_COLLECTION_T* spectrum_collection);

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
BOOLEAN_T get_spectrum_collection_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< The spectrum collection
  int first_scan,      ///< The first scan of the spectrum to retrieve 
  SPECTRUM_T* spectrum ///< The (empty) allocated SPECTRUM_T object
  );

/**
 * Adds a spectrum to the spectrum_collection.
 */
void add_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,
  SPECTRUM_T* spectrum); 

/**
 * Removes a spectrum from the spectrum_collection.
 */
void remove_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,
  SPECTRUM_T* spectrum); 

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/******************************************************************************/

/**
 * Instantiates a new spectrum_iterator object from spectrum_collection.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
SPECTRUM_ITERATOR_T* new_spectrum_iterator(
    SPECTRUM_COLLECTION_T* spectrum_collection);        

/**
 * Frees an allocated spectrum_iterator object.
 */
void free_spectrum_iterator(SPECTRUM_ITERATOR_T* spectrum_iterator);

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T spectrum_iterator_has_next(SPECTRUM_ITERATOR_T* spectrum_iterator);

/**
 * The basic iterator function next.
 */
SPECTRUM_T* spectrum_iterator_next(SPECTRUM_ITERATOR_T* spectrum_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
