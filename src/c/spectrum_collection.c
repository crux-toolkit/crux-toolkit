/*****************************************************************************
 * \file spectrum_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2003
 * DESCRIPTION: code to support working with collection of multiple spectra
 * REVISION: $Revision: 1.2 $
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

/**
 * \struct spectrum_collection 
 */
struct spectrum_collection {
  SPECTRUM_T* spectra;  ///< The spectrum peaks
  int  num_spectra;     ///< The number of peaks
  char*   filename;     ///< Optional filename
  char* comment;        ///< A comment (e.g. the spectrum_collection header lines)
  BOOLEAN_T* is_parsed; ///< Have we parsed all the spectra from the file?
};    

/**
 * \struct spectrum_iterator
 */
struct spectrum_iterator {
  SPECTRUM_COLLECTION_T* spectra;  ///< The spectrum_collection whose spectra to iterate over
  int  spectrum_idx;     ///< The index of the current spectrum
}

