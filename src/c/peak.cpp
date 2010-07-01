/*************************************************************************//**
 * \file peak.cpp
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 6/14/04
 * \brief Object for representing one peak in a spectrum.
 ****************************************************************************/
#include "peak.h"
#include "utils.h"
#include <string.h>
#include "objects.h"
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"

#include <algorithm>

using namespace std;

/**
 * \struct peak
 * \brief A spectrum peak.
 */
struct peak {
  FLOAT_T intensity;  ///< The intensity of the peak.
  FLOAT_T intensity_rank;  ///< The rank intensity of the peak.
  FLOAT_T location;   ///< The location of the peak.
};    

/**
 * \returns A PEAK_T object, 
 * heap allocated, must be freed using free_peak method
 */
PEAK_T* new_peak (
  FLOAT_T intensity, ///< intensity for the new peak -in 
  FLOAT_T location ///< location for the new peak -in
  )
{
  PEAK_T* fresh_peak =(PEAK_T*)mymalloc(sizeof(PEAK_T));
  fresh_peak->intensity = intensity;
  fresh_peak->intensity_rank = 0.0;
  fresh_peak->location = location;
  return fresh_peak;
}
/**
 * frees A PEAK_T object
 */
void free_peak (
  PEAK_T* garbage_peak ///< the peak to free -in
  )
{
  free(garbage_peak);
}

/**
 * \returns the intensity of PEAK_T object
 */
FLOAT_T get_peak_intensity(
 PEAK_T* working_peak///< return the intensity of this peak -in
) {
  return working_peak->intensity;
}

/**
 * \returns the intensity of PEAK_T object
 */
FLOAT_T get_peak_intensity_rank(
  PEAK_T* working_peak///< return the intensity of this peak -in
)
{
  return working_peak->intensity_rank;
}

/**
 * \returns the location of PEAK_T object
 */
FLOAT_T get_peak_location(
  PEAK_T* working_peak///< return the location of this peak -in 
  )
{
  return working_peak->location;
}


/**
 * sets the intensity of PEAK_T object
 */
void set_peak_intensity(
  PEAK_T* working_peak, ///< set the intensity of this peak -out
  FLOAT_T intensity///< the intensity -in
 )
{
  working_peak->intensity = intensity;
}

/**
 * sets the intensity rank of the PEAK_T object
 */
void set_peak_intensity_rank(
  PEAK_T* working_peak, ///< set the intensity rank of this peak -out
  FLOAT_T intensity_rank  ///< the intensity rank -in
 )
{
  working_peak->intensity_rank = intensity_rank;
}

/**
 * sets the location of PEAK_T object
 */
void set_peak_location(
  PEAK_T* working_peak, ///<set the location of this peak -out
  FLOAT_T location ///< the location -in
  )
{
  working_peak->location = location;
}

/**
 * prints the intensity and location of PEAK_T object to stdout
 */
void print_peak(
  PEAK_T* working_peak ///< print this peak -in
  )
{
  printf("%.1f %.1f\n", 
         working_peak->location,
         working_peak->intensity);
}

/**
 * \returns A heap allocated PEAK_T object array
 */
PEAK_T* allocate_peak_array(
  int num_peaks///< number of peaks to allocate -in
  )
{
  PEAK_T* peak_array;
  peak_array = (PEAK_T*)mycalloc(1,sizeof(PEAK_T) * num_peaks);
  return peak_array;
}


/**
 * \returns A vector of allocated PEAK_T objects
 */
vector<PEAK_T*> allocate_peak_vector(
  unsigned int num_peaks///< number of peaks to allocate -in
  )
{
  vector<PEAK_T*> ans;

  for (unsigned int idx=0;idx < num_peaks;idx++) {
    ans.push_back(new_peak(0,0));
  }

  return ans;
}

/**
 * \frees A PEAK_T object array
 */
void free_peak_array(
  PEAK_T* garbage_peak ///<the peak array to free -in
  ) 
{
  free(garbage_peak);
}

/**
 * \frees a PEAK_T object vector
 */
void free_peak_vector(
  vector<PEAK_T*>& peaks ///<the peak vector to free -in
  )
{
  for (unsigned int idx=0;idx < peaks.size();idx++) {
    free_peak(peaks[idx]);
  }
  peaks.clear();
}

/**
 *\returns a pointer to the peak in the peak_array
 */
PEAK_T* find_peak(
  PEAK_T* peak_array,///< peak_array to search -in
  int index ///< the index of peak to fine -in
  )
{
  return &peak_array[index];
}


/***********************************************
 * Sort peaks
 * also functions for lib. function sort,
 * although maybe used for other purposes
 ************************************************/
/**
 * Written for the use of lib. function, sort()
 * compare the intensity of peaks
 *\returns true if peak_1 is larger, false otherwise.
 */
bool compare_peaks_by_intensity(
  const PEAK_T* peak_one,
  const PEAK_T* peak_two
)
{
  return (peak_one->intensity > peak_two->intensity);
}

/**
 * Written for the use of lib. function, sort()
 * compare the mz(location) of peaks
 *\returns true if peak_2 is larger, false otherwise
 */
bool compare_peaks_by_mz(
  const PEAK_T* peak_one, ///< peak one to compare -in
  const PEAK_T* peak_two  ///< peak two to compare -in
  ) {

  return peak_one->location < peak_two->location;
}

/**
 * sort peaks by their intensity or location
 * use the lib. function, qsort()
 * sorts intensity in descending order
 */
void sort_peaks(
  vector<PEAK_T*>& peak_array, ///< peak array to sort -in/out
  PEAK_SORT_TYPE_T sort_type ///< the sort type(location or intensity)
  )
{
  if(sort_type == _PEAK_INTENSITY){
    // sort the peaks by intensity
    sort(peak_array.begin(), peak_array.end(), compare_peaks_by_intensity);
  }
  else if(sort_type == _PEAK_LOCATION){
    // sort the peaks by location
    sort(peak_array.begin(), peak_array.end(), compare_peaks_by_mz);
  }
  else{
    carp(CARP_ERROR, "no matching peak sort type");
  }
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
