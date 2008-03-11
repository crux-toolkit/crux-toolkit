/**
 * \file hit.c
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 2008 March 11
 * DESCRIPTION: \brief Object for collecting the evidence for a particular 
 *                     protein hit.
 * REVISION: $Revision: 1.1 $
 ****************************************************************************/
#ifndef HIT_H
#define HIT_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "objects.h"
#include "parameter.h"
#include "spectrum_collection.h"
// #include "hit_collection.h"
#include "protein.h"

/**
 * \returns a new memory allocated hit
 */
HIT_T* new_hit(void);

/**
 * free the memory allocated hit
 */
void free_hit(
  HIT_T* hit ///< the hit to free -in
  );

/**
 * print the information of the hit
 */
void print_hit(
  HIT_T* hit, ///< the hit to print -in  
  FILE* file  ///< output stream -out
);

/**
 * Increments the pointer count to the hit object
 */
void increment_hit_pointer_count(
  HIT_T* hit ///< the hit to work -in  
  );


/**
 * \returns the hit score in the hit object
 */
double get_hit_score(
  HIT_T* hit          ///< the hit to work -in  
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
