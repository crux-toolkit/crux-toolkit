/**
 * \file spectrum.h 
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 6/15/06
 * VERSION: $Revision: 1.1 $
 * DESCRIPTION: Object for representing one spectrum.
 *****************************************************************************/
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include "utils.h"
#include "peak.h"
#include <stdio.h>

/**
 * \typedef The spectrum type definition.
 */
typedef struct spectrum SPECTRUM_T;

/**
 * \returns A spectrum object.
 */
SPECTRUM_T* allocate_spectrum(void);

/**
 * \struct spectrum The spectrum struct.
 */
struct spectrum* new_spectrum(
  char*   filename,
  float mz,
  int    assumed_z,
  float  assumed_m,
  float  min_peak_mz,
  float  max_peak_mz,
  int    num_peaks,
  double total_energy,
  double* peak_mz,
  double* peak_intensity,
	int    first_scan,
	int    last_scan,
  char*   comment);

void free_spectrum
(SPECTRUM_T* spectrum);

int get_num_peaks
(SPECTRUM_T* spectrum);

float get_total_energy
(SPECTRUM_T* spectrum);

float get_mass
(SPECTRUM_T* spectrum);

float get_precursor_mz
(SPECTRUM_T* spectrum);

int get_charge
(SPECTRUM_T* spectrum);

int get_first_scan
(SPECTRUM_T* spectrum);

void parse_spectrum
(char*       filename,
 SPECTRUM_T* spectrum);

void parse_spectrum_fh
(FILE*       infile,
 SPECTRUM_T* spectrum);

void print_spectrum
(FILE*       file,
 char*       label,
 SPECTRUM_T* s);

void copy_spectrum
(SPECTRUM_T* src,
 SPECTRUM_T* dest);

#endif
