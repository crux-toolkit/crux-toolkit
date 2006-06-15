/**
 * \file peak.h
 * $Revision: 1.1 $
 * \brief Object for representing one peak in a spectrum.
 *
 * A peak is primarily identified via its intensity (height) and location
 * (position on the m/z axis).
 *
 * A peptide of length n is indexed from the N-terminal to C-terminal
 * end, with indices of 0 to n+1.  A theoretical peak, corresponding
 * to a sub-peptide, is identified by the indices of the cleavages
 * that created that sub-peptide within the original peptide.  Thus,
 * an N-terminal ion (b ion) that is created by cleaving at position
 * b4 would have an n_cleavage of 0 and a c_cleavage of 4.
 */
#ifndef PEAK_H
#define PEAK_H
#include <stdio.h>

/**
 * \typedef peak
 */
typedef struct peak {
  float intensity;
  float location;
} PEAK_T;

/**
 * \returns A PEAK_T object
 */
PEAK_T* new_peak
(float intensity,
float location
);


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
