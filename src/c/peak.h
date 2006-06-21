/**
 * \file peak.h
 * $Revision: 1.3 $
 * \brief Object for representing one peak in a spectrum.
 *
 * A peak is primarily identified via its intensity (height) and location
 * (position on the m/z axis).
 *
 */
#ifndef PEAK_H
#define PEAK_H
#include <stdio.h>

/**
 * \typedef PEAK_T 
 */
typedef struct peak {
  float intensity;
  float location;
} PEAK_T;

/**
 * \returns A PEAK_T object
 */
PEAK_T* new_peak (
  float intensity,
  float location
);

/**
 * \frees A PEAK_T object
 */
void free_peak (PEAK_T* garbage_peak);

/**
 * \returns the intensity of PEAK_T object
 */
float peak_intensity(PEAK_T* working_peak);

/**
 * \returns the location of PEAK_T object
 */
float peak_location(PEAK_T* working_peak);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
