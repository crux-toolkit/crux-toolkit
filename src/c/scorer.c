/*****************************************************************************
 * \file scorer.c
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * DESCRIPTION: object to score spectrum vs. spectrum or spectrum vs. ion_series
 * REVISION: $Revision: 1.9 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "spectrum.h"
#include "scorer.h"
#include "parameter.h"


//the bin width
#define bin_width_mono 1.0005079
#define bin_width_average 1.0011413

/**
 * \struct scorer
 * \brief An object to score spectrum v. spectrum or spectrum v. ion_series
 */
struct scorer {
  SCORER_TYPE_T type; ///< The type of scorer
  float sp_beta; ///< used for Sp: the beta variable 
  //float sp_gamma; ///< used for Sp: the gamma variable 
  //float sp_min_mz; ///< used for Sp: the min mz for the intensity sum array
  float sp_max_mz; ///< used for Sp: the max mz for the intensity array
  //float sp_array_resolution; ///< used for Sp: the resolution for intensity sum array indecies
  //float sp_sum_resolution; ///<  used for Sp: the resolution to which interval of peaks to add of intensity sum
  //float sp_equalize_resolution; ///<  used for Sp: the resolution to which interval of peaks should be equalized to the highest peak
  float* intensity_array; ///< used for Sp: the intensity array, which can be indexed using the m/z
  float max_intensity; ///< the max intensity in the intensity array
};

/**
 *\returns An (empty) scorer object.
 */
SCORER_T* allocate_scorer(void){
  SCORER_T* scorer = (SCORER_T*)mycalloc(1, sizeof(SCORER_T));
  return scorer;
}

/**
 * If not planning to use the default values, must parse the parameter file before.
 * Instantiates a new scorer object from a filename. 
 * \returns a new scorer object
 */
SCORER_T* new_scorer(
  SCORER_TYPE_T type ///< the type of scorer -in
  )
{
  //get allocated scorer
  SCORER_T* scorer = allocate_scorer();
  
  //set score type
  scorer->type = type;
  
  //set fields needed for each score type
  if(type == SP){
    scorer->sp_beta = get_double_parameter("beta", 0.075);
    //scorer->sp_gamma = get_double_parameter("gama", 0.15);
    //scorer->sp_min_mz = get_double_parameter("min_mz", 0); 
    scorer->sp_max_mz = get_double_parameter("max-mz", 4000);
    //scorer->sp_array_resolution = get_double_parameter("sp_array_resolution", 0.5);
    //scorer->sp_sum_resolution = get_double_parameter("sp_sum_resolution", 1);
    //scorer->sp_equalize_resolution = get_double_parameter("sp_equalize_resolution", 1);

    //allocate the intensity array
    scorer->intensity_array = (float*)mycalloc(scorer->sp_max_mz, sizeof(float));
    scorer->max_intensity = 0;
  }

  printf("%s\n",get_string_parameter("dummy"));
  
  return scorer;
}

/**
 * Frees an allocated scorer object.
 */
void free_scorer(
  SCORER_T* scorer ///< the ion collection to free - in
  )
{
  //score type SP?
  if(scorer->type == SP){
    free(scorer->intensity_array);
  }

  free(scorer);
}


/**
 * normalize array so that maximum peak equals threshold
 */
void nomalize_intensity_array(
  float* intensity_array, ///< the array to normalize -in/out
  int array_size, ///< size of array -in
  float max_intensity, ///< the maximum intensity in array -in
  float threshold   ///< the threshold to which the peaks should be normalized -in
  )
{
  int mz_idx = 0;

  //return if max_intensity is 0
  if(max_intensity < 0.00001){
    return;
  }

  //normalize all peaks
  for(; mz_idx < array_size; ++mz_idx){
    intensity_array[mz_idx] = intensity_array[mz_idx] * threshold / max_intensity;

    //DEBUG
    carp(CARP_INFO, "norm data[%d] = %.2f",mz_idx, intensity_array[mz_idx]); 
  }
}

/**
 * smooth all peaks in intensity array
 * Replaces the original array with the newly smooothed array
 */
void smooth_peaks(
  SCORER_T* scorer        ///< the scorer object -in/out
  )
{
  int idx = 2;
  float* array = scorer->intensity_array;

  //create a new array, which will replace the original intensity array
  float* new_array = (float*)mycalloc(scorer->sp_max_mz, sizeof(float));

  
  switch (scorer->type){
    case SP:
      //iterate over all peaks
      for(; idx < (int)scorer->sp_max_mz-2; ++idx){
        //smooooth
        new_array[idx] = (array[idx-2]+4*array[idx-1]+6*array[idx]+4*array[idx+1]+array[idx+2])/16;

        //DEBUG
        carp(CARP_INFO, "smooth data[%d] = %f",idx, new_array[idx]); 
      }
      break;
      
  case XCORR:
    break;
    
  case DOTP:
    break;
  }
  free(scorer->intensity_array);
  scorer->intensity_array = new_array;
}

/**
 * get the mean of intensity in array within +/- 50 mz of the working peak
 * \returns the mean +/- 50mz region
 */
float get_mean_from_array(
  float* original_array, ///< the array to normalize -in
  int array_size, ///< the size of array -in
  int peak_idx,  ///< the peak indx in array -in
  int* peak_count ///< pointer to peak count, store peak count here -out
  )
{
  float total_intensity = 0;
  int start_idx = peak_idx - 50;
  int end_idx = peak_idx + 50;

  //set upper bound
  if(peak_idx + 50 >= array_size){
    end_idx = array_size-1;
  }
  //set start index
  if(peak_idx - 50 <= 0){
    start_idx = 0;
  }
  
  //sum up the intensities
  for(; start_idx <= end_idx; ++start_idx){
    //printf("%.2f\n", original_array[start_idx]);

    ++*peak_count;
    total_intensity += original_array[start_idx];
  }
  
  //BUG! it should divide by 101 but Sequest uses 100
  return (total_intensity / (*peak_count-1));
}

/**
 * get the stdev of intensity in array within +/- 50 mz of the working peak
 * \returns the stdev +/- 50mz region
 */
float get_stdev_from_array(
  float* original_array, ///< the array to normalize -in
  int array_size, ///< the size of array -in
  int peak_idx,  ///< the peak indx in array -ina
  float mean,     ///< the mean in the +/- 50 interval -in
  int peak_count ///<  peak count -in
  )
{
  float variance = 0;
  int start_idx = peak_idx - 50;
  int end_idx = peak_idx + 50;
  float dev = 0;

  //set upper bound
  if(peak_idx + 50 >= array_size){
    end_idx = array_size-1;
  }
  //set start index
  if(peak_idx - 50 <= 0){
    start_idx = 0;
  }
  
  //sum up the intensities
  for(; start_idx <= end_idx; ++start_idx){
    //sum up all deviations squared
    dev = original_array[start_idx] - mean;
    variance += (dev*dev);
  }
  
  //return the stdev
  return sqrt(variance/peak_count);
}

/***
 * zero and extract peaks
 * extract peaks that are larger than mean + #step*stdev into new array
 * zero out the peaks that have been extracted
 * yes, the facter that a peak has removed will effect the fallowing peaks
 */
void zero_peak_mean_stdev(
  SCORER_T* scorer,        ///< the scorer object -in/out
  float* original_array, ///< the array to normalize -in/out
  float* new_array, ///< the array to normalize -in/out                          
  int step                ///< is this 1 or 2 step -in
  )
{
  int peak_count = 0;
  int idx = 0;
  int array_size = (int)scorer->sp_max_mz;
  float mean = 0;
  float stdev = 0;

  //iterate over all peaks
  for(; idx < array_size; ++idx){
    peak_count = 0;
    //get mean
    mean = get_mean_from_array(original_array, array_size, idx, &peak_count);
    //get stdev
    stdev = get_stdev_from_array(original_array, array_size, idx, mean, peak_count);
    
    //DEBUG
    carp(CARP_INFO, "zero idx: %d mean: %.8f, stdev: %.8f", idx, mean, stdev);
    
    //iterate over all positions and extract peaks
    if(original_array[idx] > (mean + step*stdev)){
      //set new array with modified intensity
      new_array[idx] = original_array[idx] - (mean - stdev);
      
      //DEBUG
      carp(CARP_INFO, "extract peak: %.2f at idx: %d", original_array[idx], idx);
      
      //for step 1,
      if(step == 1){
        //zero out original peak
        original_array[idx] = 0;
      }
    }
  }
}

/**
 *  zero and extract peaks
 * extract peaks that are larger than mean + #step*stdev into new array
 * zero out the peaks that have been extracted
 * repeat twice, than replace old array with extracted peak array
 */
void zero_peaks(
  SCORER_T* scorer   ///< the scorer object -in/out
  )
{
  //create a new array, which will replace the original intensity array
  float* new_array = (float*)mycalloc(scorer->sp_max_mz, sizeof(float));
  
  //step 1,
  zero_peak_mean_stdev(scorer, scorer->intensity_array, new_array, 1);
  //step 2,
  zero_peak_mean_stdev(scorer, scorer->intensity_array, new_array, 2);
  
  //replace intensity_array with new intensity array
  free(scorer->intensity_array);
  scorer->intensity_array = new_array;
}

/**
 * keep only the peaks up to top rank peaks remove other peaks.
 * do second normalization on the top peaks back to max 100 intensity
 * replace old array with normalized top peak array
 */
void extract_peaks(
  SCORER_T* scorer,        ///< the scorer object -in/out
  int top_rank  ///< keep the top ranking peaks -in
  )
{
  //create a new array, which will replace the original intensity array
  float* temp_array = (float*)mycalloc((int)scorer->sp_max_mz, sizeof(float));
  float* original_array = scorer->intensity_array;
  int idx = 0;
  int temp_idx = 0;
  float cut_off = 0;
  float max_intensity;

  //copy all peaks to temp_array
  for(; idx < (int)scorer->sp_max_mz; ++idx){
    if(scorer->intensity_array[idx] > 0){
      temp_array[temp_idx] = original_array[idx];

      //DEBUG print all temp array values
      carp(CARP_INFO, "before sort data[%d]=%.3f",temp_idx, temp_array[temp_idx]);
      
      ++temp_idx;
    }
  }
  
  //if there's over top_rank peaks, keep only top_rank peaks
  //quick sort
  quicksort(temp_array, temp_idx);
  
  //set max and cut_off
  max_intensity = temp_array[0];
  cut_off = temp_array[top_rank-1] + 0.0001;
  
  //remove peaks bellow cut_off 
  //also, normalize peaks to max_intensity to 100
  for(idx = 0; idx < (int)scorer->sp_max_mz; ++idx){
    //DEBUG print all temp array values
    carp(CARP_INFO, "scorted data[%d]=%.3f",idx, temp_array[idx]);

    if(original_array[idx] > 0){
      //is it bellow cut off?
      if(original_array[idx] < cut_off){
        //remove peak
        original_array[idx] = 0;
      }
      else{
        //nomalize peak to max 100
        original_array[idx] = original_array[idx] / max_intensity * 100;
      }
    }
  }
  
  free(temp_array);
  return;
}
 
/**
 * equalize peaks +/- 1 mz to largest peak
 * start from left to right
 * yes, the former peaks equalization will affect the later peaks
 */
void equalize_peaks(
  SCORER_T* scorer ///< the scorer object -in/out
  )
{
  int idx;

  //equalize peaks to it's greatest intensity +/- window
  for(idx = 1; idx < (int)scorer->sp_max_mz-1; ++idx){
    if(scorer->intensity_array[idx] > 0){
    
      if(scorer->intensity_array[idx] < scorer->intensity_array[idx-1]){
        scorer->intensity_array[idx] = scorer->intensity_array[idx-1];
      }
      
      if(scorer->intensity_array[idx] < scorer->intensity_array[idx+1]){
        scorer->intensity_array[idx] = scorer->intensity_array[idx+1];
      }

      //DEBUG print all temp array values
      carp(CARP_INFO, "final data[%d]=%.3f",idx, scorer->intensity_array[idx]);
    }
  }
}
    
/**
 * create the intensity array
 * SCORER must have been created for SP type
 * \returns TRUE if successful, else FLASE
 */
BOOLEAN_T create_intensity_array(
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  SCORER_T* scorer        ///< the scorer object -in/out
  )
{
  PEAK_T* peak = NULL;
  PEAK_ITERATOR_T* peak_iterator = NULL;
  float max_intensity = 0;
  int mz = 0;
  float intensity = 0;
  float experimental_mass_cut_off = get_spectrum_precursor_mz(spectrum)*get_int_parameter("charge",2) + 50;
  int precursor_mz = (int)(get_spectrum_precursor_mz(spectrum)+0.5);

  //FIXME, later be able pick between average and mono
  float bin_width = bin_width_mono;

  //DEBUG
  carp(CARP_INFO, "precursor_mz: %d", precursor_mz);
  
  //if score type equals SP
  if(scorer->type != SP){
    carp(CARP_ERROR, "incorrect scorer type, only use this method for SP scorers");
    return FALSE;
  }
  
  //check bin width type
  //DO this at some time!!

  //create a peak iterator
  peak_iterator = new_peak_iterator(spectrum);
  
  //while there are more peaks to iterate over..
  while(peak_iterator_has_next(peak_iterator)){
    peak = peak_iterator_next(peak_iterator);
    mz = (int)(get_peak_location(peak)/bin_width + 0.5);
    
    //skip all peaks larger than experimental mass
    if(mz > experimental_mass_cut_off){
      continue;
    }
    
    //skip all peaks within precursor ion mz +/- 15
    if(mz <= precursor_mz + 15 && mz >= precursor_mz - 15){
      continue;
    }

    //get intensity
    intensity = sqrt(get_peak_intensity(peak));
    
    //set intensity in array with correct mz, only if max peak in the bin
    if(scorer->intensity_array[mz] < intensity){
      scorer->intensity_array[mz] = intensity;
      
      //check if max_intensity
      //TEST ME, THIS COULD CHAGE!!
      if(intensity > max_intensity){
        max_intensity = intensity;
      }
    }
  }

  //set max_intensity
  scorer->max_intensity = max_intensity;
  
  //DEBUG!!
  carp(CARP_INFO, "exp_cut_off: %.2f max_intensity: %.2f", experimental_mass_cut_off, (max_intensity*max_intensity));

  
  //normalize intensity
  nomalize_intensity_array(scorer->intensity_array, scorer->sp_max_mz, scorer->max_intensity, 100);
  
  //smooth peaks
  smooth_peaks(scorer);

  //zero peaks
  zero_peaks(scorer);
  
  //keep top 200
  extract_peaks(scorer, 200);

  //equalize peaks
  equalize_peaks(scorer);
  
  
  free_peak_iterator(peak_iterator);
  return TRUE;
}

/**
 * calculates all the necessay values for Sp score, related to the specfic ion_type
 * adds to intensity_sum and repeat_count
 *\returns the number of matches found from the predicted ions
 */
int calculate_ion_type_sp(
  SCORER_T* scorer,        ///< the scorer object -in                          
  ION_SERIES_T* ion_series, ///< the ion series to score against the spectrum -in
  float* intensity_sum,     ///< the total intensity sum of all matches so far -out
  ION_TYPE_T ion_type,      ///< the ion type to check -in
  int* repeat_count         ///< the repeated count of ions (ex. consecutive b ions) -out
  )
{
  int cleavage_idx = 0;
  ION_T* ion = NULL;
  float one_intensity = 0;
  int ion_match = 0;
  int ion_charge = 0;

  int* before_cleavage = (int*)mymalloc(get_ion_series_charge(ion_series)*sizeof(int));
  int cleavage_array_idx = 0;

  //initialize before cleavage indecies
  for(; cleavage_array_idx < get_ion_series_charge(ion_series); ++cleavage_array_idx){
    before_cleavage[cleavage_array_idx] = -1;
  }

  //create ion constraint
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint(get_ion_constraint_mass_type(get_ion_series_ion_constraint(ion_series)), get_ion_series_charge(ion_series), ion_type, FALSE);
  
  //create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = new_ion_filtered_iterator(ion_series, ion_constraint);
  
  //while there are ion's in ion iterator, add matched observed peak intensity
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ion = ion_filtered_iterator_next(ion_iterator);
    
    //get the intensity matching to ion's m/z
    one_intensity = scorer->intensity_array[(int)(get_ion_mass_z(ion)/bin_width_mono + 0.5)];
    
    //if there is a match in the observed spectrum
    if(one_intensity > 0){
      
      //DEBUG
      carp(CARP_INFO, "matched ion: %.2f ion intensity: %.2f", get_ion_mass_z(ion), one_intensity);

      ++ion_match;
      *intensity_sum = *intensity_sum + one_intensity;
      
      //get ion charge
      ion_charge = get_ion_charge(ion) - 1;
      
      //check if repeated ion b1, b2, ...
      if((cleavage_idx = get_ion_cleavage_idx(ion)) == before_cleavage[ion_charge] + 1){
        ++*repeat_count;
      }
      
      //reset the previous cleavage index 
      before_cleavage[ion_charge] = cleavage_idx;
    }
  }
  
  
  //free ion iterator, ion_constraint
  free(before_cleavage);
  free_ion_constraint(ion_constraint);
  free_ion_filtered_iterator(ion_iterator);

  return ion_match;
}


/**
 * given a spectrum and ion series calculates the Sp score
 *\returns the sp score 
 */
float gen_score_sp(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  ION_SERIES_T* ion_series ///< the ion series to score against the spectrum -in
  )
{
  float final_score = 0;
  float intensity_sum = 0;
  int ion_match = 0;
  int repeat_count = 0;
    
  //get gama and beta
  float beta = scorer->sp_beta;

  //create intensity array
  if(!create_intensity_array(spectrum, scorer)){
    carp(CARP_ERROR, "failed to produce Sp");
    free(spectrum);
    free(ion_series);
    free(scorer);
    exit(1);
  }
  
  //calculate the B_ION and Y_ION portions of the Sp score
  ion_match = calculate_ion_type_sp(scorer, ion_series, &intensity_sum, B_ION, &repeat_count) +
    calculate_ion_type_sp(scorer, ion_series, &intensity_sum, Y_ION, &repeat_count);
  
  //// DEBUG!!!!
  carp(CARP_INFO, "# repeat count: %d ion_match count: %d, total_ion count: %d sum: %.2f", 
       repeat_count, ion_match, get_ion_series_num_ions(ion_series),
       intensity_sum);

  //calculate Sp score.
  if(ion_match != 0){
    final_score = 
      (intensity_sum * ion_match) * (1+ (repeat_count * beta)) / get_ion_series_num_ions(ion_series);
  }
  
  //return score
  return final_score;
}

/**
 * Score a spectrum vs. an ion series
 */
float score_spectrum_v_ion_series(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  ION_SERIES_T* ion_series ///< the ion series to score against the spectrum -in
  )
{
  float final_score = 0;

  //if score type equals SP
  if(scorer->type == SP){
    final_score = gen_score_sp(scorer, spectrum, ion_series);
  }
  //FIXME, later add different score types...
  else{
    carp(CARP_ERROR, "no scoring method availiable for the scorers' score type");
  }
  
  return final_score;
}

/**
 * Score a spectrum vs. another spectrum
 */
float score_spectrum_v_spectrum(
  SCORER_T* scorer,           ///< the scorer object
  SPECTRUM_T* first_spectrum, ///< the first spectrum to score 
  SPECTRUM_T* second_spectrum ///<  the second spectrum to score
);



/*******************************
 * get, set methods for scorer
 *******************************/

/**
 *\returns the score type of the scorer
 */
SCORER_TYPE_T get_scorer_type(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->type;
}

/**
 *\sets the scorer type
 */
void set_scorer_type(
  SCORER_T* scorer, ///< the scorer object -out                     
  SCORER_TYPE_T type ///< The type of scorer -in
)
{
  scorer->type = type;
}

/**
 *\returns the beta value of the scorer
 */
float get_scorer_sp_beta(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_beta;
}

/**
 *\sets the scorer beta value
 */
void set_scorer_sp_beta(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_beta ///< used for Sp: the beta variable -in
  )
{
  scorer->sp_beta = sp_beta;
}

/**
 *\returns the gamma value of the scorer
 */
/*
float get_scorer_sp_gamma(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_gamma;
}
*/

/**
 *\sets the scorer gamma value
 */
/*
void set_scorer_sp_gamma(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_gamma ///< used for Sp: the gamma variable -in
  )
{
  scorer->sp_gamma = sp_gamma;
}
*/

/**
 *\returns the min_mz value of the scorer
 */
/*
float get_scorer_sp_min_mz(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_min_mz;
}
*/

/**
 *\sets the scorer min_mz value
 */
/*
void set_scorer_sp_min_mz(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_min_mz ///< used for Sp: the min_mz variable -in
  )
{
  scorer->sp_min_mz = sp_min_mz;
}
*/

/**
 *\returns the max_mz value of the scorer
 */
float get_scorer_sp_max_mz(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_max_mz;
}

/**
 *\sets the scorer max_mz value
 */
void set_scorer_sp_max_mz(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_max_mz ///< used for Sp: the max_mz variable -in
  )
{
  scorer->sp_max_mz = sp_max_mz;
}

/**
 *\returns the array_resolution value of the scorer
 */
/*
float get_scorer_sp_array_resolution(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_array_resolution;
}
*/

/**
 *\sets the scorer array_resolution value
 */
/*
void set_scorer_sp_array_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_array_resolution ///< used for Sp: the array_resolution variable -in
  )
{
  scorer->sp_array_resolution = sp_array_resolution;
}
*/

/**
 *\returns the sum_resolution value of the scorer
 */
/*
float get_scorer_sp_sum_resolution(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_sum_resolution;
}
*/

/**
 *\sets the scorer sum_resolution value
 */
/*
void set_scorer_sp_sum_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_sum_resolution ///< used for Sp: the sum_resolution variable -in
  )
{
  scorer->sp_sum_resolution = sp_sum_resolution;
}
*/

/**
 *\returns the equalize_resolution value of the scorer
 */
/*
float get_scorer_sp_equalize_resolution(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_equalize_resolution;
}
*/

/**
 *\sets the scorer equalize_resolution value
 */
/*
void set_scorer_sp_equalize_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_equalize_resolution ///< used for Sp: the equalize_resolution variable -in
  )
{
  scorer->sp_equalize_resolution = sp_equalize_resolution;
}
*/


/**************************
 *OLD code
 *
 **************************/

/**
 * calculates all the necessay values for Sp score, related to the specfic ion_type
 * adds to intensity_sum and repeat_count
 *\returns the number of matches found from the predicted ions
 */
/*
int calculate_ion_type_sp(
  SCORER_T* scorer,        ///< the scorer object -in                          
  SPECTRUM_T* spectrum,     ///< the spectrum to score -in
  ION_SERIES_T* ion_series, ///< the ion series to score against the spectrum -in
  float* intensity_sum,     ///< the total intensity sum of all matches so far -out
  ION_TYPE_T ion_type,      ///< the ion type to check -in
  int* repeat_count         ///< the repeated count of ions (ex. consecutive b ions) -out
  )
{
  int before_cleavage_idx = -1;
  int cleavage_idx = 0;
  ION_T* ion = NULL;
  int one_intensity = 0;
  int ion_match = 0;

  //create ion constraint
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint(get_ion_constraint_mass_type(get_ion_series_ion_constraint(ion_series)), 1, ion_type, FALSE);
  
  //create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = new_ion_filtered_iterator(ion_series, ion_constraint);
  
  //while there are ion's in ion iterator, add matched observed peak intensity
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ion = ion_filtered_iterator_next(ion_iterator);
    one_intensity = get_nearby_intensity_sum(scorer, spectrum, get_ion_mass_z(ion));
    
    //if there is a match in the observed spectrum
    if(one_intensity != 0){
      ++ion_match;
      *intensity_sum += one_intensity;
      
      //check if repeated ion b1, b2, ...
      if((cleavage_idx = get_ion_cleavage_idx(ion)) == before_cleavage_idx + 1){
        ++*repeat_count;
      }
      
      //reset the previous cleavage index 
      before_cleavage_idx = cleavage_idx;
    }
  }
    
  //free ion iterator, ion_constraint
  free_ion_constraint(ion_constraint);
  free_ion_filtered_iterator(ion_iterator);

  return ion_match;
}
*/

/**
 * given a spectrum and ions calculates the Sp score
 *\returns the sp score 
 */
/*
// old SP implementation..
float gen_score_sp(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  ION_SERIES_T* ion_series ///< the ion series to score against the spectrum -in
  )
{
  float final_score = 0;
  float intensity_sum = 0;
  int ion_match = 0;
  int repeat_count = 0;
  int amino_count = 0;
  SPECTRUM_T* processed_spectrum = NULL;

  //get gama and beta
  float gamma = scorer->sp_gamma; 
  float beta = scorer->sp_beta;
  
  //get a processed spectrum
  processed_spectrum = process_spectrum(spectrum, SP);

  //calculate the B_ION and Y_ION portions of the Sp score
  ion_match = calculate_ion_type_sp(scorer, processed_spectrum, ion_series, &intensity_sum, B_ION, &repeat_count) +
    calculate_ion_type_sp(scorer, processed_spectrum, ion_series, &intensity_sum, Y_ION, &repeat_count);
  
  //normalize the intensity_sum to 100
  intensity_sum = intensity_sum * 100.0 / get_spectrum_max_intensity(processed_spectrum);

  //// DEBUG!!!!
  //carp(CARP_INFO, "# max intensity: %.2f ion_match count: %d, total_ion count: %d", 
  //     get_spectrum_max_intensity(processed_spectrum), ion_match, get_ion_series_num_ions(ion_series));

  //calculate Sp score.
  if(ion_match != 0){
    final_score = 
      intensity_sum * ion_match * (1+repeat_count * beta) * (1+amino_count * gamma) / get_ion_series_num_ions(ion_series);
  }
  
  //free the processed_spectrum
  free_spectrum(processed_spectrum);

  //return score
  return final_score;
}
*/

/**
 * given a spectrum and ion series calculates the Sp score
 *\returns the sp score 
 */
/*
float gen_score_sp(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  ION_SERIES_T* ion_series ///< the ion series to score against the spectrum -in
  )
{
  float final_score = 0;
  float intensity_sum = 0;
  int ion_match = 0;
  int repeat_count = 0;
  SPECTRUM_T* processed_spectrum = NULL;
  
  //get gama and beta
  float beta = scorer->sp_beta;

  //get a processed spectrum
  processed_spectrum = process_spectrum(spectrum, SP);

  //fill in the intensity array, that is normalized to 100 in scorer
  if(!fill_intensity_array(processed_spectrum, scorer)){
    carp(CARP_ERROR, "failed to produce Sp");
    free(spectrum);
    free(ion_series);
    free(scorer);
    exit(1);
  }
  
  //calculate the B_ION and Y_ION portions of the Sp score
  ion_match = calculate_ion_type_sp(scorer, ion_series, &intensity_sum, B_ION, &repeat_count) +
    calculate_ion_type_sp(scorer, ion_series, &intensity_sum, Y_ION, &repeat_count);
  
  //// DEBUG!!!!
  //  carp(CARP_INFO, "# max intensity: %.2f ion_match count: %d, total_ion count: %d sum: %.2f", 
  //   get_spectrum_max_intensity(processed_spectrum), ion_match, get_ion_series_num_ions(ion_series),
  //   intensity_sum);

  //calculate Sp score.
  if(ion_match != 0){
    final_score = 
      (intensity_sum * ion_match) * (1+ (repeat_count * beta))) / get_ion_series_num_ions(ion_series);
    //(get_ion_series_charge(ion_series) 
  }
  
//free the processed_spectrum
  free_spectrum(processed_spectrum);

  //return score
  return final_score;
}
*/
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
