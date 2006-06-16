/**
 * \file spectrum.h 
 * $Revision: 1.2 $
 * \brief Object for representing one spectrum.
 *****************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "utils.h"
#include "peak.h"
#include <stdio.h>

/**
 * \typedef SPECTRUM_T 
 */
typedef struct spectrum SPECTRUM_T;

/**
 * \typedef SPECTRUM_TYPE_T The spectrum type (MS, MS-MS, MS-MS-MS)
 */
typedef enum type { ms1, ms2, ms3 } SPECTRUM_TYPE_T;

/**
 * \returns An (empty) spectrum object.
 */
SPECTRUM_T* allocate_spectrum(void);

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */
SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan
	int               last_scan,          ///< The number of the last scan
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum.
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z,         ///< The possible charge states of this spectrum
  char*             filename);          ///< Optional filename

/**
 * Frees an allocated spectrum object.
 */
void free_spectrum (SPECTRUM_T* spectrum);

/**
 * Prints a spectrum object to file.
 */
void print_spectrum(SPECTRUM_T* spectrum, FILE* file);

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(SPECTRUM_T* spectrum);

/**
 * Copies spectrum object src to dest.
 */
void copy_spectrum(
  SPECTRUM_T* src,
  SPECTRUM_T* dest);

/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 */
int parse_spectrum_file(
  SPECTRUM_T* spectrum,
  FILE* file);

/**
 * \returns TRUE if success. FALSE is failure.
 */
int parse_spectrum(
  SPECTRUM_T* spectrum,
  char*      filename);

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * Each module should also provide new_<object>, free_<object> and print_<object> functions.
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 * \returns The intensity of the peak with the maximum intensity.
 */
float get_spectrum_base_peak_intensity(SPECTRUM_T* spectrum);

/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
float get_spectrum_mass(SPECTRUM_T* spectrum, int charge);

/**
 * \returns The mass of the neutral precursor ion, according to the formula 
 * mass = m/z * charge - mass_H * charge
 */
float get_spectrum_neutral_mass(SPECTRUM_T* spectrum, int charge);

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
float get_spectrum_singly_charged_mass(SPECTRUM_T* spectrum);

/**
 * Adds a peak to the spectrum, modifying num_peaks, and max_ and min_peak_mz, if necessary.
 */
float add_peak(SPECTRUM_T* spectrum, PEAK_T* peak);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
