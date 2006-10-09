/**
 * \file scorer.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * $Revision: 1.2 $
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
 */
SCORER_T* new_scorer(
  SCORER_TYPE_T type,
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
  SCORER_T* scorer,        ///< the scorer object
  SPECTRUM_T* spectrum,    ///< the spectrum to score
  ION_SERIES_T* scorer     ///< the ion series to score against the spectrum
);

/**
 * Score a spectrum vs. another spectrum
 */
float score_spectrum_v_spectrum(
  SCORER_T* scorer,           ///< the scorer object
  SPECTRUM_T* first_spectrum, ///< the first spectrum to score 
  SPECTRUM_T* second_spectrum ///<  the second spectrum to score
);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
