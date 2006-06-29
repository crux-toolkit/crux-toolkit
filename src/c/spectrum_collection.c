/*****************************************************************************
 * \file spectrum_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2003
 * DESCRIPTION: code to support working with collection of multiple spectra
 * REVISION: $Revision: 1.5 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "spectrum.h"
#include "peak.h"
#include "utils.h"

#define MAX_SPECTRA 40000

/**
 * \struct spectrum_collection 
 */
struct spectrum_collection {
  SPECTRUM_T* spectra[MAX_SPECTRA];  ///< The spectrum peaks
  int  num_spectra;     ///< The number of spectra
  char*   filename;     ///< Optional filename
  char* comment;        ///< The spectrum_collection header lines
  BOOLEAN_T is_parsed; ///< Have we parsed all the spectra from the file?
};    

/**
 * \struct spectrum_iterator
 */
struct spectrum_iterator {
  SPECTRUM_T* spectrum; ///< The spectrum whose peaks to iterate over. 
  int  peak_index;        ///< The index of the current peak
};


/**
 * \returns An (empty) spectrum_collection object.
 */
SPECTRUM_COLLECTION_T* allocate_spectrum_collection(void){
  SPECTRUM_COLLECTION_T* collection =
    (SPECTRUM_COLLECTION_T*)mycalloc(1,sizeof(SPECTRUM_COLLECTION_T));
  collection->is_parsed = FALSE;
  return collection;
}

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does *NOT* parse all spectra in the file. 
 * This will be done lazily depending on the subsequent method
 * calls (parse_spectrum_collection get_spectrum_collection_spectrum).
 * \returns  SPECTRUM_COLLECTION_T
 */
SPECTRUM_COLLECTION_T* new_spectrum_collection(char* filename)///< The spectrum collection filename.
{
  SPECTRUM_COLLECTION_T* spectrum_collection =  allocate_spectrum_collection();

  if(!access(filename, F_OK)){
    fprintf(stderr,"File %s could not be opened",filename);
    exit(1);
  } //FIXME check if file is empty
  
  set_spectrum_collection_filename(spectrum_collection, filename);
  
  return spectrum_collection;
}

/**
 * Frees an allocated spectrum_collection object.
 */
void free_spectrum_collection(SPECTRUM_COLLECTION_T* spectrum_collection){
  int spectrum_index = 0;
  for(;spectrum_index < spectrum_collection->num_spectra; ++spectrum_index){
    free_spectrum(spectrum_collection->spectra[spectrum_index]);
  }
  free(spectrum_collection->filename);
  free(spectrum_collection->commet);
  free(spectrum_collection);
}


/**  //////TESTME////
 * \sets the filename of the ms2 file the spectra was parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_collection_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_collection_filename(SPECTRUM_COLLECTION_T* spectrum_collection, char* filename){
  int filename_length = strlen(filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  spectrum_collection->filename =
    strncpy(copy_filename,filename,filename_length);  
}


/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan. Use binary search
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
BOOLEAN_T get_spectrum_collection_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< The spectrum collection
  int first_scan,      ///< The first scan of the spectrum to retrieve 
  SPECTRUM_T* spectrum ///< The (empty) allocated SPECTRUM_T object
  )
{
  FILE* file;
  long low_index;
  long high_index;

  if ((file = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"File %s could not be opened",filename);
    return (FALSE);
  }
  
  low_index = ftell(file);
  


  return TRUE;
}

long binary_search_spectrum(long low_index, long high_index);



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
