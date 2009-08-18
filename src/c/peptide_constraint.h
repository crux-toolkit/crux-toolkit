/**
 * \file peptide_constraint.h 
 * $Revision: 1.5 $
 * \brief Object for holding the peptide constraint information.
 */
#ifndef PEPTIDE_CONSTRAINT_H 
#define PEPTIDE_CONSTRAINT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "crux-utils.h"
#include "objects.h"
#include "mass.h"
#include "peptide.h"
#include "protein.h"
#include "carp.h"


/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* allocate_peptide_constraint(void);

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint(
  ENZYME_T enzyme, ///< the enzyme to use for digestion
  DIGEST_T digest, ///< the degree of digestion
//  PEPTIDE_TYPE_T peptide_type, ///< the peptide_type -in
  FLOAT_T min_mass, ///< the minimum mass -in
  FLOAT_T max_mass, ///< the maximum mass -in
  int min_length, ///< the minimum length of peptide -in
  int max_length, ///< the maximum lenth of peptide -in
  int num_mis_cleavage, ///< The maximum mis cleavage of the peptide -in
  MASS_TYPE_T mass_type  ///< isotopic mass type (AVERAGE, MONO) -in
  );

/**
 * \brief Create a new peptide constraint and populate its values
 * based on those in parameter.c 
 * \returns A newly allocated peptide constraint.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint_from_parameters();

/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T peptide_constraint_is_satisfied(
   PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraints to enforce -in
   PEPTIDE_T* peptide ///< the query peptide -in
   );

/**
 * Copies an allocated peptide_constraint object.
 */
PEPTIDE_CONSTRAINT_T* copy_peptide_constraint_ptr(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< object to copy -in 
  );

/**
 * Frees an allocated peptide_constraint object.
 */
void free_peptide_constraint(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< object to free -in 
  );


/**peptide_constraint**/

/**
 * sets the peptide type of the peptide_constraint
 */
/*
void set_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide constraint - in
  );
*/
/**
 * \returns the peptide type of the peptide_constraint
 */
/*
PEPTIDE_TYPE_T get_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );
*/
void set_peptide_constraint_enzyme(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  ENZYME_T enzyme
);

ENZYME_T get_peptide_constraint_enzyme(
  PEPTIDE_CONSTRAINT_T* peptide_constraint///< the peptide constraint to set -out
);

void set_peptide_constraint_digest(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  DIGEST_T digest
);

DIGEST_T get_peptide_constraint_digest(
  PEPTIDE_CONSTRAINT_T* peptide_constraint///< the peptide constraint to set -out
);

/**
 * sets the min mass of the peptide_constraint
 */
void set_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide constraint to set -out
  FLOAT_T min_mass  ///< the min mass of the peptide constraint - in
  );

/**
 * \returns the min mass of the peptide_constraint
 */
FLOAT_T get_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the max mass of the peptide_constraint
 */
void set_peptide_constraint_max_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  FLOAT_T max_mass  ///< the max mass of the peptide constraint - in
  );

/**
 * \returns the max mass of the peptide_constraint
 */
FLOAT_T get_peptide_constraint_max_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the min length of the peptide_constraint
 */
void set_peptide_constraint_min_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int min_length  ///< the min length of the peptide constraint - in
  );

/**
 * \returns the min length of the peptide_constraint
 */
int get_peptide_constraint_min_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the max length of the peptide_constraint
 */
void set_peptide_constraint_max_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int max_length  ///< the max length of the peptide constraint - in
  );

/**
 * \returns the max length of the peptide_constraint
 */
int get_peptide_constraint_max_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );


/**
 * sets the num_mis_cleavage of the peptide_constraint
 */
void set_peptide_constraint_num_mis_cleavage(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int num_mis_cleavage ///< The maximum mis cleavage of the peptide -in
  );

/**
 * \returns the num_mis_cleavage of the peptide_constraint
 */
int get_peptide_constraint_num_mis_cleavage(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the mass type of the peptide_constraint
 */
void set_peptide_constraint_mass_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  MASS_TYPE_T mass_type ///< the peptide_type for the constraint -in
  );

/**
 * \returns the mass type of the mass_constraint
 */
MASS_TYPE_T get_peptide_constraint_mass_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

#ifdef __cplusplus
}
#endif

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
