/**
 * \file scorer.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * $Revision: 1.5 $
 * DESCRIPTION: object to score spectrum vs. spectrum or spectrum vs. scorer
 *****************************************************************************/
#ifndef SCORER_H 
#define SCORER_H

#include <stdio.h>
#include "objects.h"
#include "peptide.h"
#include "ion.h"
#include "scorer.h"

/**
 * \returns An (empty) scorer object.
 */
SCORER_T* allocate_scorer(void);

/**
 * Instantiates a new scorer object from a filename. 
 * \returns a new scorer object
 */
SCORER_T* new_scorer(
  SCORER_TYPE_T type ///< the type of scorer -in
  );

/**
 * Frees an allocated scorer object.
 */
void free_scorer(
  SCORER_T* scorer ///< the ion collection to free - in
);

/**
 * Score a spectrum vs. an ion series
 */
float score_spectrum_v_ion_series(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum,    ///< the spectrum to score -in
  ION_SERIES_T* ion_series ///< the ion series to score against the spectrum -in
);

/**
 * Score a spectrum vs. another spectrum
 */
float score_spectrum_v_spectrum(
  SCORER_T* scorer,           ///< the scorer object -in
  SPECTRUM_T* first_spectrum, ///< the first spectrum to score -in
  SPECTRUM_T* second_spectrum ///<  the second spectrum to score -in
);


/*******************************
 * get, set methods for scorer
 *******************************/

/**
 *\returns the score type of the scorer
 */
SCORER_TYPE_T get_scorer_type(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer type
 */
void set_scorer_type(
  SCORER_T* scorer, ///< the scorer object -out                     
  SCORER_TYPE_T type ///< The type of scorer -in
  );

/**
 *\returns the beta value of the scorer
 */
float get_scorer_sp_beta(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer beta value
 */
void set_scorer_sp_beta(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_beta ///< used for Sp: the beta variable -in
  );

/**
 *\returns the gamma value of the scorer
 */
float get_scorer_sp_gamma(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer gamma value
 */
void set_scorer_sp_gamma(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_gamma ///< used for Sp: the gamma variable -in
  );


/**
 *\returns the min_mz value of the scorer
 */
float get_scorer_sp_min_mz(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer min_mz value
 */
void set_scorer_sp_min_mz(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_min_mz ///< used for Sp: the min_mz variable -in
  );


/**
 *\returns the max_mz value of the scorer
 */
float get_scorer_sp_max_mz(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer max_mz value
 */
void set_scorer_sp_max_mz(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_max_mz ///< used for Sp: the max_mz variable -in
  );


/**
 *\returns the sp_array_resolution value of the scorer
 */
float get_scorer_sp_array_resolution(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer sp_array_resolution value
 */
void set_scorer_sp_array_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_array_resolution ///< used for Sp: the sp_array_resolution variable -in
  );

/**
 *\returns the sp_sum_resolution value of the scorer
 */
float get_scorer_sp_sum_resolution(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer sp_sum_resolution value
 */
void set_scorer_sp_sum_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_sum_resolution ///< used for Sp: the sp_sum_resolution variable -in
  );

/**
 *\returns the equalize_resolution value of the scorer
 */
float get_scorer_sp_equalize_resolution(
  SCORER_T* scorer ///< the scorer object -in
  );

/**
 *sets the scorer equalize_resolution value
 */
void set_scorer_sp_equalize_resolution(
  SCORER_T* scorer, ///< the scorer object -out                     
  float sp_equalize_resolution ///< used for Sp: the equalize_resolution variable -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
