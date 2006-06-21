/*****************************************************************************
 * \file spectrum.c
 * AUTHOR: Tobias P. Mann
 * CREATE DATE: 19 Sept 2003
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.3 $
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
 * \struct spectrum 
 */
struct spectrum {
  int               first_scan;   ///< The number of the first scan
  int               last_scan;    ///< The number of the last scan
  int               id;           ///< A unique identifier
  SPECTRUM_TYPE_T   spectrum_type;///< The type of spectrum.
  float             precursor_mz; ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z;   ///< The possible charge states of this spectrum
  PEAK_T*           peaks;        ///< The spectrum peaks
  float             min_peak_mz;  ///< The minimum m/z of all peaks
  float             max_peak_mz;  ///< The maximum m/z of all peaks
  int               num_peaks;    ///< The number of peaks
  double            total_energy; ///< The sum of intensities in all peaks
  char*             filename;     ///< Optional filename
};    

/**
 * \returns An (empty) spectrum object.
 */
SPECTRUM_T* allocate_spectrum(){
  SPECTRUM_T* fresh_spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));
  return fresh_spectrum;
}

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */
SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan
  int               last_scan,          ///< The number of the last scan
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum.
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z,         ///< The possible charge states of this spectrum
  char*             filename)          ///< Optional filename
{
  SPECTRUM_T* fresh_spectrum = allocate_spectrum();
  fresh_spectrum->first_scan = first_scan;
  fresh_spectrum->last_scan = last_scan;
  fresh_spectrum->spectrum_type = spectrum_type;
  fresh_spectrum->precursor_mz = precursor_mz;
  fresh_spectrum->possible_z = possible_z;
  fresh_spectrum->filename = filename; //optional, needs some work done on this.is it the file it reads?
}

/**
 * Frees an allocated spectrum object.
 */
void free_spectrum (SPECTRUM_T* spectrum){
  free(spectrum->possible_z);
  free_peak(spectrum->peaks);
  free(spectrum->filename);
  free(spectrum);
}


/**
 * Prints a spectrum object to file.
 */
void print_spectrum(SPECTRUM_T* spectrum, FILE* file){

}

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(SPECTRUM_T* spectrum){


}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T
 */
void copy_spectrum(
  SPECTRUM_T* src,
  SPECTRUM_T* dest)
{
  //copy each varible
  //copy each peak

}


/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_spectrum_file(
  SPECTRUM_T* spectrum,
  FILE* file)
{


}

/**
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_spectrum(
  SPECTRUM_T* spectrum,
  char*      filename)
{



}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

