/*****************************************************************************
 * \file peak.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 6/14/04
 * VERSION: $Revision: 1.7 $
 * DESCRIPTION: Object for representing one peak in a spectrum.
 *****************************************************************************/
#include "peak.h"
#include "utils.h"
#include <string.h>

/**
 * \struct peak
 * \brief A spectrum peak.
 */
struct peak {
  float intensity;  ///< The intensity of the peak.
  float location;   ///< The location of the peak.
};    

/**
 * \returns A PEAK_T object, 
 * heap allocated, must be freed using free_peak method
 */
PEAK_T* new_peak (float intensity, float location){
  PEAK_T* fresh_peak =(PEAK_T*)mymalloc(sizeof(PEAK_T));
  fresh_peak->intensity = intensity;
  fresh_peak->location = location;
  return fresh_peak;
}
/**
 * frees A PEAK_T object
 */
void free_peak (PEAK_T* garbage_peak){
  free(garbage_peak);
}

/**
 * \returns the intensity of PEAK_T object
 */
float get_peak_intensity(PEAK_T* working_peak){
  return working_peak->intensity;
}

/**
 * \returns the location of PEAK_T object
 */
float get_peak_location(PEAK_T* working_peak){
  return working_peak->location;
}


/**
 * sets the intensity of PEAK_T object
 */
void set_peak_intensity(PEAK_T* working_peak, float intensity){
  working_peak->intensity = intensity;
}

/**
 * sets the location of PEAK_T object
 */
void set_peak_location(PEAK_T* working_peak, float location){
  working_peak->location = location;
}

/**
 * prints the intensity and location of PEAK_T object to stdout
 */
void print_peak(PEAK_T* working_peak){
  printf("%.1f %.1f\n", 
         working_peak->location,
         working_peak->intensity);
}

/**
 * \returns A heap allocated PEAK_T object array
 */
PEAK_T* allocate_peak_array(int num_peaks){
  PEAK_T* peak_array;
  peak_array = (PEAK_T*)mycalloc(1,sizeof(PEAK_T) * num_peaks);
  return peak_array;
}

/**
 * \frees A PEAK_T object array
 */
void free_peak_array(PEAK_T* garbage_peak){
  free(garbage_peak);
}

/**
 *\returns a pointer to the peak in the peak_array
 */
PEAK_T* find_peak(PEAK_T* peak_array, int index){
  return &peak_array[index];
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
