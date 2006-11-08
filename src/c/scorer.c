/*****************************************************************************
 * \file scorer.c
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * DESCRIPTION: object to score spectrum vs. spectrum or spectrum vs. ion_series
 * REVISION: $Revision: 1.7 $
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
    scorer->sp_max_mz = get_double_parameter("max_mz", 4000);
    //scorer->sp_array_resolution = get_double_parameter("sp_array_resolution", 0.5);
    //scorer->sp_sum_resolution = get_double_parameter("sp_sum_resolution", 1);
    //scorer->sp_equalize_resolution = get_double_parameter("sp_equalize_resolution", 1);

    //allocate the intensity array
    scorer->intensity_array = (float*)mycalloc(scorer->sp_max_mz, sizeof(float));
  }
  
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
 * fills in the intensity array
 * SCORER must have been created for SP type
 * \returns TRUE if successful, else FLASE
 */
BOOLEAN_T fill_intensity_array(
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  SCORER_T* scorer        ///< the scorer object -in
  )
{
  PEAK_T* peak = NULL;
  PEAK_ITERATOR_T* peak_iterator = NULL;
  float max_intensity = 0;
  int mz = 0;
  float intensity = 0;

  //if score type equals SP
  if(scorer->type != SP){
    carp(CARP_ERROR, "incorrect scorer type, only use this method for SP scorers");
    return FALSE;    
  }
  
  //get max_intensity to normalize peaks to 100
  max_intensity = get_spectrum_max_intensity(spectrum);

  //create a peak iterator
  peak_iterator = new_peak_iterator(spectrum);
  
  //while there are more peaks to iterate over..
  while(peak_iterator_has_next(peak_iterator)){
    peak = peak_iterator_next(peak_iterator);
    mz = (int)get_peak_location(peak);
    
    //normalize peak intensity
    intensity = get_peak_intensity(peak) * 100.0 / max_intensity;

    //set intensity in array with correct mz
    scorer->intensity_array[mz] = intensity;
  }
  
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
  int before_cleavage_idx = -1;
  int cleavage_idx = 0;
  ION_T* ion = NULL;
  int one_intensity = 0;
  int ion_match = 0;

  //create ion constraint
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint(get_ion_constraint_mass_type(get_ion_series_ion_constraint(ion_series)), get_ion_series_charge(ion_series), ion_type, FALSE);
  
  //create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = new_ion_filtered_iterator(ion_series, ion_constraint);
  
  //while there are ion's in ion iterator, add matched observed peak intensity
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ion = ion_filtered_iterator_next(ion_iterator);
    
    //DEBUG
    /*
    if((int)get_ion_mass_z(ion) == 728){
      continue;
    }
    */

    //print_ion(ion, stdout);
    //FIXME, maybe don't add 0.5 to m/z..
    //get the intensity matching to ion's m/z
    one_intensity = scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))] +
      scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+1] + 
      scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-1] +
      scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+2] + 
      scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-2];
     

    //DEBUG
    float temp = scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))];
    
    /*
    int temp = 0;
    int n = 0;
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))]){
      ++n;
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+1]){
      ++n;
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-1]){
      ++n;
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+2]){
      ++n;
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-2]){
      ++n;
    }
    if(n>0){
      one_intensity = one_intensity/n;
    }
    */
    
   
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+1]){
      temp =  scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+1];
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-1]){
      temp =  scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-1];
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+2]){
      temp =  scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+2];
    }
    if(temp < scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-2]){
      temp =  scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-2];
    }
    one_intensity = temp;
    

    
    //if there is a match in the observed spectrum
    if(one_intensity != 0){

      //DEBUG
      //printf("matched ion: %.2f\n", get_ion_mass_z(ion));

      ++ion_match;
      *intensity_sum = *intensity_sum + one_intensity;// +
        //scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))+1] + 
        //scorer->intensity_array[(int)(0.5 + get_ion_mass_z(ion))-1];
            
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
  SPECTRUM_T* processed_spectrum = NULL;
  
  //get gama and beta
  //float beta = scorer->sp_beta;

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
      (intensity_sum * ion_match) /* (1+ (repeat_count * beta)))*/  / get_ion_series_num_ions(ion_series);
    //(get_ion_series_charge(ion_series) 
  }
  
  //free the processed_spectrum
  free_spectrum(processed_spectrum);

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


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
