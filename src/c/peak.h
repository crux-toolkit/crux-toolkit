/**
 * \file peak.h
 * $Revision: 1.14 $
 * \brief Object for representing one peak in a spectrum.
 *
 * A peak is primarily identified via its intensity (height) and location
 * (position on the m/z axis).
 *
 */
#ifndef PEAK_H
#define PEAK_H
#include <stdio.h>
#include <stdlib.h>
#include "objects.h"
#include "utils.h"

#include <vector>

/**
 * \returns A PEAK_T object
 */
PEAK_T* new_peak (
  FLOAT_T intensity, ///< intensity for the new peak -in 
  FLOAT_T location ///< location for the new peak -in
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
FLOAT_T get_peak_intensity(
  PEAK_T* working_peak ///< return the intensity of this peak -in
  );

/**
 * sets the intensity rank of PEAK_T object
 */
FLOAT_T get_peak_intensity_rank(
  PEAK_T* working_peak ///< get the intensity rank of this peak -in
  );

/**
 * \returns the location of PEAK_T object
 */
FLOAT_T get_peak_location(
  PEAK_T* working_peak ///< return the location of this peak -in 
  );

/**
 * sets the intensity of PEAK_T object
 */
void set_peak_intensity(
  PEAK_T* working_peak, ///< set the intensity of this peak -mod
  FLOAT_T intensity ///< the intensity -in
  );

/**
 * sets the intensity rank of PEAK_T object
 */
void set_peak_intensity_rank(
  PEAK_T* working_peak, ///< set the intensity of this peak -mod
  FLOAT_T intensity_rank ///< the intensity -in
  );

/**
 * sets the location of PEAK_T object
 */
void set_peak_location(
  PEAK_T* working_peak, ///<set the location of this peak -out
  FLOAT_T location ///< the location -in
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
 * \returns A vector of allocated PEAK_T objects
 */
std::vector<PEAK_T*> allocate_peak_vector(
  unsigned int num_peaks///< number of peaks to allocate -in
  );

/**
 * \frees A PEAK_T object array
 */
void free_peak_array(
  PEAK_T* garbage_peak ///<the peak array to free -in
  );

/**
 * \frees a PEAK_T object vector
 */
void free_peak_vector(
  std::vector<PEAK_T*>& peaks ///<the peak array to free -in
  );

/**
 *\returns a pointer to the peak in the peak_array
 */
PEAK_T* find_peak(
  PEAK_T* peak_array,///< peak_array to search -in
  int index ///< the index of peak to fine -in
  );

/***********************************************
 * Sort peaks
 * also functions for lib. function sort(),
 * although maybe used for other purposes
 ************************************************/

/**
 * Written for the use of lib. function, sort()
 * compare the intensity of peaks
 *\returns true if peak_1 is larger, false otherwise.
 */
bool compare_peaks_by_intensity(
  const PEAK_T* peak_one, ///< peak one to compare -in
  const PEAK_T* peak_two  ///< peak two to compare -in
  );

/**
 * Written for the use of lib. function, sort()
 * compare the mz(location) of peaks
 *\returns true if peak_2 is larger, false otherwise
 */
bool compare_peaks_by_mz(
  const PEAK_T* peak_one, ///< peak one to compare -in
  const PEAK_T* peak_two  ///< peak two to compare -in
  );

/**
 * sort peaks by their intensity or location
 * use the lib. function, sort()
 */
void sort_peaks(
  std::vector<PEAK_T*>& peak_array, ///< peak array to sort -in/out
  PEAK_SORT_TYPE_T sort_type ///< the sort type(location or intensity)
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
