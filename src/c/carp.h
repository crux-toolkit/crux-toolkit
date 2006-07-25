/**
 * \file carp.h 
 * $Revision: 1.1 $
 * \brief Provides methods for logging error messages, and setting verbosity level.
 *****************************************************************************/
#ifndef CARP_H 
#define CARP_H

#define CARP_FATAL 0
#define CARP_ERROR 10
#define CARP_WARNING 20
#define CARP_INFO 30
#define CARP_DEBUG 40
#define CARP_DETAILED_DEBUG 50
#define CARP_MAX 60

#include <stdio.h>
#include "utils.h"

/**
 * Sets the global verbosity level. Expects an integer greater than or equal
 * to CARP_FATAL and less than CARP_MAX
 */
void set_verbosity_level(int verbosity);

/**
 * \returns the current verbosity level.
 */
int get_verbosity_level(void);

/**
 * \returns True if error message has successfully been output. False if not.
 * The term 'carp' is used because 'log' is already used by the math library. 
 */
BOOLEAN_T carp(
    int verbosity,
    char* format, 
    ...);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
