/*****************************************************************************
 * \file scorer.c
 * AUTHOR: Chris Park
 * CREATE DATE: 9 Oct 2006
 * DESCRIPTION: object to score spectrum vs. spectrum or spectrum vs. ion_series
 * REVISION: $Revision: 1.2 $
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

/**
 * \struct scorer
 * \brief An object to score a series of ions, and organize them!
 * For which additional data structures will be created as needed 
 */
struct scorer {
  SCORER_TYPE_T type; ///< The type of scorer
};




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
