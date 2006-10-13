/*****************************************************************************
 * \file scorer.c
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * DESCRIPTION: object to score spectrum vs. spectrum or spectrum vs. ion_series
 * REVISION: $Revision: 1.4 $
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
  float sp_gamma; ///< used for Sp: the gamma variable 
  float sp_min_mz; ///< used for Sp: the min mz for the intensity sum array
  float sp_max_mz; ///< used for Sp: the max mz for the intensity sum array
  float sp_array_resolution; ///< used for Sp: the resolution for intensity sum array indecies
  float sp_sum_resolution; ///<  used for Sp: the resolution to which interval of peaks to add of intensity sum
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
    scorer->sp_gamma = get_double_parameter("gama", 0.15);
    scorer->sp_min_mz = get_double_parameter("min_mz", 0); 
    scorer->sp_max_mz = get_double_parameter("max_mz", 4000);
    scorer->sp_array_resolution = get_double_parameter("sp_array_resolution", 0.5);
    scorer->sp_sum_resolution = get_double_parameter("sp_sum_resolution", 1);
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
  free(scorer);
}

/**
 * calculates all the necessay values for Sp score, related to the specfic ion_type
 * adds to intensity_sum and repeat_count
 *\returns the number of matches found from the predicted ions
 */
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

  //create ion constraint
  ION_CONSTRAINT_T* ion_constraint = 
    new_ion_constraint(get_ion_constraint_mass_type(get_ion_series_constraint(ion_series)), 1, ion_type, FALSE);
  
  //create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = new_ion_filtered_iterator(ion_series, ion_constraint);
  
  //while there are ion's in ion iterator, add matched observed peak intensity
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ion = ion_filtered_iterator_next(ion_iterator);
    one_intensity = get_nearby_intensity_sum(scorer, specturm, get_ion_mz(ion));
    
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
    
  //free ion iterator
  free_ion_filtered_iterator(ion_iterator);
  
  return ion_match;
}

/**
 * given a spectrum and ions calculates the Sp score
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
  float gamma = scorer->sp_gamma; 
  float beta = scorer->sp_beta;
  
  //calculate the B_ION and Y_ION portions of the Sp score
  int ion_match = calculate_ion_type_sp(scorer, spectrum, ion_series, &intensity_sum, B_ION, &repeat_count) +
    calculate_ion_type_sp(scorer, spectrum, ion_series, &intensity_sum, Y_ION, &repeat_count);
  
  //calculate Sp score.
  final_score = 
    intensity_sum * ion_match * (1+repeat_count * beta) * (1+amino_count * gama) / get_ion_series_num_ions(ion_series);
  
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
  if(scorer->scorer_type == SP){
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
float get_scorer_sp_gamma(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_gamma;
}

/**
 *\sets the scorer gamma value
 */
void set_scorer_sp_gamma(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_gamma ///< used for Sp: the gamma variable -in
  )
{
  scorer->sp_gamma = sp_gamma;
}


/**
 *\returns the min_mz value of the scorer
 */
float get_scorer_sp_min_mz(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_min_mz;
}

/**
 *\sets the scorer min_mz value
 */
void set_scorer_sp_min_mz(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_min_mz ///< used for Sp: the min_mz variable -in
  )
{
  scorer->sp_min_mz = sp_min_mz;
}

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
float get_scorer_sp_array_resolution(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_array_resolution;
}

/**
 *\sets the scorer array_resolution value
 */
void set_scorer_sp_array_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_array_resolution ///< used for Sp: the array_resolution variable -in
  )
{
  scorer->sp_array_resolution = sp_array_resolution;
}

/**
 *\returns the sum_resolution value of the scorer
 */
float get_scorer_sp_sum_resolution(
  SCORER_T* scorer ///< the scorer object -in
  )
{
  return scorer->sp_sum_resolution;
}

/**
 *\sets the scorer sum_resolution value
 */
void set_scorer_sp_sum_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_sum_resolution ///< used for Sp: the sum_resolution variable -in
  )
{
  scorer->sp_sum_resolution = sp_sum_resolution;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
