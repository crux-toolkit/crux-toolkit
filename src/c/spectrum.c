/*****************************************************************************
 * \file spectrum.c
 * AUTHOR: Tobias P. Mann
 * CREATE DATE: 19 Sept 2003
 * DESCRIPTION: code to support working with spectra
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


