/**
 * \file peak.h
 * $Revision: 1.10 $
 * \brief Object for representing one peak in a spectrum.
 *
 * A peak is primarily identified via its intensity (height) and location
 * (position on the m/z axis).
 *
 */
#ifndef PEAK_H
#define PEAK_H
#include <stdio.h>
#include "objects.h"

/**
 * \returns A PEAK_T object
 */
PEAK_T* new_peak (
  float intensity, ///< intensity for the new peak -in 
  float location ///< location for the new peak -in
  );

/**
 * \frees A PEAK_T object
 */
void free_peak (
  PEAK_T* garbage_peak ///< the peak to free -in
  );

/**
 * \returns the intensity of PEAK_T object
 */
float get_peak_intensity(
  PEAK_T* working_peak///< return the intensity of this peak -in
  );

/**
 * \returns the location of PEAK_T object
 */
float get_peak_location(
  PEAK_T* working_peak///< return the location of this peak -in 
  );


/**
 * sets the intensity of PEAK_T object
 */
void set_peak_intensity(
  PEAK_T* working_peak, ///<set the intensity of this peak -out
  float intensity///< the intensity -in
  );

/**
 * sets the location of PEAK_T object
 */
void set_peak_location(
  PEAK_T* working_peak, ///<set the location of this peak -out
  float location ///< the location -in
  );

/**
 * \prints the intensity and location of PEAK_T object to stdout
 */
void print_peak(
  PEAK_T* working_peak ///< print this peak -in
  );

/**
 * \returns A heap allocated PEAK_T object array
 */
PEAK_T* allocate_peak_array(
  int num_peaks///< number of peaks to allocate -in
  );


/**
 * \frees A PEAK_T object array
 */
void free_peak_array(
  PEAK_T* garbage_peak ///<the peak array to free -in
  ); 

/**
 *\returns a pointer to the peak in the peak_array
 */
PEAK_T* find_peak(
  PEAK_T* peak_array,///< peak_array to search -in
  int index ///< the index of peak to fine -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
