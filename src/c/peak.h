/**
 * \file peak.h
 * $Revision: 1.8 $
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
 * A peak in a spectrum
 */
typedef struct peak PEAK_T;

/**
 * \returns A PEAK_T object
 */
PEAK_T* new_peak (
  float intensity,
  float location);

/**
 * \frees A PEAK_T object
 */
void free_peak (PEAK_T* garbage_peak);

/**
 * \returns the intensity of PEAK_T object
 */
float get_peak_intensity(PEAK_T* working_peak);

/**
 * \returns the location of PEAK_T object
 */
float get_peak_location(PEAK_T* working_peak);

/**
 * sets the intensity of PEAK_T object
 */
void set_peak_intensity(PEAK_T* working_peak, float intensity);

/**
 * sets the location of PEAK_T object
 */
void set_peak_location(PEAK_T* working_peak, float location);

/**
 * \prints the intensity and location of PEAK_T object to stdout
 */
void print_peak(PEAK_T* working_peak);


/**
 * \returns A heap allocated PEAK_T object array
 */
PEAK_T* allocate_peak_array(int num_peaks);


/**
 * \frees A PEAK_T object array
 */
void free_peak_array(PEAK_T* garbage_peak);

/**
 *\returns a pointer to the peak in the peak_array
 */
PEAK_T* find_peak(PEAK_T* peak_array, int index);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
