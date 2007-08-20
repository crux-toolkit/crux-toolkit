/**
 * \file scorer.h 
 * \brief object to score spectrum vs. spectrum or spectrum vs. scorer
 */

/*
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * $Revision: 1.14 $
 *****************************************************************************/
#ifndef SCORER_H 
#define SCORER_H

#include <stdio.h>
#include <dirent.h>
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
 * Frees the single_ion_constraints array
 */
void free_single_ion_constraints(
  ION_CONSTRAINT_T** ion_constraints
);

/**
 * Creates the an array of ion constraints for GMTK models.
 * TODO do we need one for paired and single? Do we want an iterator?
 */
ION_CONSTRAINT_T** single_ion_constraints(
    void
);

/**
 * Create ion files (for GMTK) in the output directory for all psms for a
 * single spectrum. Appends to existing files.
 * \returns TRUE for success
 */
BOOLEAN_T output_psm_files(
  char* output_directory,  ///< name of directory to place the ion files
  SPECTRUM_T* spectrum,    ///< input spectrum
  char** peptides, ///< the peptide sequences
  int num_peptides ///< the number of peptides sequences
);

/**
 * Score a spectrum vs. another spectrum
 */
float score_spectrum_v_spectrum(
  SCORER_T* scorer,           ///< the scorer object -in
  SPECTRUM_T* first_spectrum, ///< the first spectrum to score -in
  SPECTRUM_T* second_spectrum ///<  the second spectrum to score -in
);

/*************************************
 * Score for LOGP_*
 ************************************/

/**
 * Compute a p-value for a given score w.r.t. an exponential with the given parameters.
 *\returns the -log(p_value) of the exponential distribution
 */
float score_logp_exp_sp(
  float sp_score, ///< The sp score for the scoring peptide -in
  float mean      ///< The overall mean of the sp scored peptides -in
  );

/**
 * Compute a p-value for a given score w.r.t. an exponential with the given parameters.
 *\returns the -log(p_value) of the exponential distribution with Bonferroni correction
 */
float score_logp_bonf_exp_sp(
  float sp_score, ///< The sp score for the scoring peptide -in
  float mean,      ///< The overall mean of the sp scored peptides -in
  int num_peptide  ///< The number of peptides scored for sp
  );

/**
 * Compute a p-value for a given score w.r.t. a Weibull with given parameters.
 *\returns the -log(p_value)
 */
float score_logp_weibull(
  float score, ///< The score for the scoring peptide -in
  float eta,  ///< The eta parameter of the Weibull
  float beta ///< The beta parameter of the Weibull
  );

/**
 * Compute a p-value for a given score w.r.t. a Weibull with given parameters.
 *\returns the -log(p_value)
 */
float score_logp_bonf_weibull(
  float score, ///< The score for the scoring peptide -in
  float eta,  ///< The eta parameter of the Weibull
  float beta, ///< The beta parameter of the Weibull
  int num_peptides ///< The number of peptides.
  );


/**
 * Compute a p-value for a given score w.r.t. an EVD with the given parameters.
 *\returns the -log(p_value) of the EVD distribution 
 */
float score_logp_evd_xcorr(
  float xcorr_score, ///< The xcorr score for the scoring peptide -in
  float mu, ///<  EVD parameter Xcorr(characteristic value of extreme value distribution) -in
  float l_value ///< EVD parameter Xcorr(decay constant of extreme value distribution) -in
  );

/**
 * Compute a p-value for a given score w.r.t. an EVD with the given parameters.
 *\returns the -log(p_value) of the EVD distribution with Bonferroni correction
 */
float score_logp_bonf_evd_xcorr(
  float xcorr_score, ///< The xcorr score for the scoring peptide -in
  float mu, ///<  EVD parameter Xcorr(characteristic value of extreme value distribution) -in
  float l_value, ///< EVD parameter Xcorr(decay constant of extreme value distribution) -in
  int num_peptide  ///< The number of peptides scored for sp -in
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

/**
 *\returns the fraction of b,y ions matched for scoring SP, the values is valid for the last ion series scored with this scorer object
 */
float get_scorer_sp_b_y_ion_match(
  SCORER_T* scorer ///< the scorer object -out
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
