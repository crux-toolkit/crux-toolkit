/**
 * \file hit.c
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 2008 March 11
 * DESCRIPTION: \brief Object for collecting the evidence for a particular 
 *                     protein hit.
 * REVISION: $Revision: 1.5 $
 ****************************************************************************/
#ifndef HIT_H
#define HIT_H

#ifdef __cplusplus
extern "C" {
#endif

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
 * \returns Increments the hit score by score
 */
void hit_increment_score(
    HIT_T* hit,
    double score);

/**
 * \returns Gets the hit protein;
 */
PROTEIN_T* get_hit_protein(
    HIT_T* hit);

/**
 * \returns Sets the hit protein
 */
void set_hit_protein(
    HIT_T* hit,
    PROTEIN_T* protein);

/**
 * performs protein level normalization of hit score
 */
void hit_recalibrate_score(
   HIT_T* hit
);


/**
 * print the information of the hit
 */
void print_hit(
  FILE* file,  ///< output stream -out
  HIT_T* hit ///< the hit to print -in  
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
#ifdef __cplusplus
}
#endif

#endif
