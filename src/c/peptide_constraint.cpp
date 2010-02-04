/*************************************************************************//**
 * \file peptide_constraint.cpp
 * $Revision: 1.10 $
 * \brief: Object for holding the peptide constraint information.
 ****************************************************************************/
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
#include "peptide_constraint.h"

/**
 * \struct peptide_constraint
 * \brief Object to represent constraints which a peptide may or may not
 *  satisfy.
 *
 * def TRYPTIC: a protein that ends with either K or R and 
 *              any other K and R in the sequence must be followed by a P
 */
struct peptide_constraint {
//  PEPTIDE_TYPE_T peptide_type;///< The type of peptides(TRYPTIC, PARTIALLY_TRYPTIC, N_TRYPTIC, C_TRYPTIC, NOT_TRYPTIC, ANY_TRYPTIC)
  ENZYME_T enzyme;
  DIGEST_T digestion;
  FLOAT_T min_mass; ///< The minimum mass of the peptide
  FLOAT_T max_mass; ///< The maximum mass of the peptide
  int min_length; ///< The minimum length of the peptide
  int max_length; ///< The maximum length of the peptide
  int num_mis_cleavage; ///< The maximum mis cleavage of the peptide
  MASS_TYPE_T mass_type; ///< isotopic mass type (AVERAGE, MONO)
  int num_pointers; ///< Number of pointers to this constraint
};



/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* allocate_peptide_constraint(void){
  PEPTIDE_CONSTRAINT_T* peptide_constraint =
    (PEPTIDE_CONSTRAINT_T*)mycalloc(1, sizeof(PEPTIDE_CONSTRAINT_T));
  peptide_constraint->num_pointers = 1;
  carp(CARP_DETAILED_DEBUG, "Allocating peptide constraint");
  return peptide_constraint;
}

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint(
//  PEPTIDE_TYPE_T peptide_type, ///< The type of peptides, is it TRYPTIC -in
  ENZYME_T enzyme, 
  DIGEST_T digest,
  FLOAT_T min_mass, ///< the minimum mass -in
  FLOAT_T max_mass, ///< the maximum mass -in
  int min_length, ///< the minimum length of peptide -in
  int max_length,  ///< the maximum lenth of peptide(max limit = 255) -in
  int num_mis_cleavage, ///< The maximum mis cleavage of the peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  )
{
  // max length must be less or equal than 255 becuase of the unsigned char limit of 255
  if(max_length > 255){
    carp(CARP_FATAL, "ERROR: cannot set max length higer than 255");
  }

  PEPTIDE_CONSTRAINT_T* peptide_constraint =
    allocate_peptide_constraint();

//  set_peptide_constraint_peptide_type(peptide_constraint, peptide_type);
  set_peptide_constraint_enzyme(peptide_constraint, enzyme);
  set_peptide_constraint_digest(peptide_constraint, digest);
  set_peptide_constraint_min_mass(peptide_constraint, min_mass);
  set_peptide_constraint_max_mass(peptide_constraint, max_mass);
  set_peptide_constraint_min_length(peptide_constraint, min_length);
  set_peptide_constraint_max_length(peptide_constraint, max_length);
  set_peptide_constraint_num_mis_cleavage(peptide_constraint, num_mis_cleavage);
  set_peptide_constraint_mass_type(peptide_constraint, mass_type);
  return peptide_constraint;
}

/**
 * \brief Create a new peptide constraint and populate its values
 * based on those in parameter.c 
 * \returns A newly allocated peptide constraint.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint_from_parameters(){
  PEPTIDE_CONSTRAINT_T* new_constraint = allocate_peptide_constraint();
  //new_constraint->peptide_type = get_peptide_type_parameter("cleavages");
  new_constraint->enzyme = get_enzyme_type_parameter("enzyme");
  new_constraint->digestion = get_digest_type_parameter("digestion");
  new_constraint->min_mass = get_double_parameter("min-mass");
  new_constraint->max_mass = get_double_parameter("max-mass");
  new_constraint->min_length = get_int_parameter("min-length");
  new_constraint->max_length = get_int_parameter("max-length");
  new_constraint->num_mis_cleavage = get_int_parameter("missed-cleavages");
  new_constraint->mass_type = get_mass_type_parameter("isotopic-mass");

  return new_constraint;
}


/**
 * Copy peptide pointer and increment pointer count
 */
PEPTIDE_CONSTRAINT_T* copy_peptide_constraint_ptr(
    PEPTIDE_CONSTRAINT_T* constraint){
  constraint->num_pointers++;
  return constraint;
}

// FIXME check the association..as long as there is one tryptic parent then true
// num_miss_cleavage is not implemented..add if needed
/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T peptide_constraint_is_satisfied(
    PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraints to enforce -in
    PEPTIDE_T* peptide ///< the query peptide -in
    )
{
  if(get_peptide_length(peptide) <= get_peptide_constraint_max_length(peptide_constraint) &&
     get_peptide_length(peptide) >= get_peptide_constraint_min_length(peptide_constraint) &&
     get_peptide_peptide_mass(peptide) <= get_peptide_constraint_max_mass(peptide_constraint) &&
     get_peptide_peptide_mass(peptide) >= get_peptide_constraint_min_mass(peptide_constraint)
     )
    {
      return TRUE;
    }
  
  return FALSE;
}

/**
 * Frees an allocated peptide_constraint object.
 */
void free_peptide_constraint(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< object to free -in 
  )
{
  peptide_constraint->num_pointers--;
  if (peptide_constraint->num_pointers == 0){
    carp(CARP_DETAILED_DEBUG, "Final free of peptide constraint");
    free(peptide_constraint);
  }
}


/**
 * sets the peptide type of the peptide_constraint
 */
/*
void set_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  PEPTIDE_TYPE_T peptide_type ///< the peptide_type for the constraint -in
  )
{
  peptide_constraint->peptide_type = peptide_type;
}
*/
/**
 * \returns the peptide type of the peptide_constraint
 */
/*
PEPTIDE_TYPE_T get_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->peptide_type;
}
*/
/**
 * Sets the enzyme used for the in silicos digestion
 * of the protein sequence into peptides.
 */
void set_peptide_constraint_enzyme(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  ENZYME_T enzyme
){
  peptide_constraint->enzyme = enzyme;
}

/**'
 * \returns The enzyme for this peptide constraint.
 */
ENZYME_T get_peptide_constraint_enzyme(
  PEPTIDE_CONSTRAINT_T* peptide_constraint///< the peptide constraint to set -out
){
  return peptide_constraint->enzyme;
}

/**
 * Sets the level of digestion for the peptide constraint.
 */
void set_peptide_constraint_digest(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  DIGEST_T digest
){
peptide_constraint->digestion = digest;
}

/**
 * \returns The level of digestion for the peptide constraint.
 */
DIGEST_T get_peptide_constraint_digest(
  PEPTIDE_CONSTRAINT_T* peptide_constraint///< the peptide constraint to set -out
){
return peptide_constraint->digestion;
}


/**
 * sets the min mass of the peptide_constraint
 */
void set_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide constraint to set -out
  FLOAT_T min_mass  ///< the min mass of the peptide constraint - in
  )
{
  peptide_constraint->min_mass = min_mass;
}

/**
 * \returns the min mass of the peptide_constraint
 */
FLOAT_T get_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->min_mass;
}


/**
 * sets the max mass of the peptide_constraint
 */
void set_peptide_constraint_max_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  FLOAT_T max_mass  ///< the max mass of the peptide constraint - in
  )
{
  peptide_constraint->max_mass = max_mass;
}

/**
 * \returns the max mass of the peptide_constraint
 */
FLOAT_T get_peptide_constraint_max_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->max_mass;
}

/**
 * sets the min length of the peptide_constraint
 */
void set_peptide_constraint_min_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int min_length  ///< the min length of the peptide constraint - in
  )
{
  peptide_constraint->min_length = min_length;
}

/**
 * \returns the min length of the peptide_constraint
 */
int get_peptide_constraint_min_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->min_length;
}

/**
 * sets the max length of the peptide_constraint
 * maximum limit 255
 */
void set_peptide_constraint_max_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int max_length  ///< the max length of the peptide constraint - in
  )
{
  // check if maximum length is with in range <= 255
  if(max_length > 255){
    carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
  }
  
  peptide_constraint->max_length = max_length;
}

/**
 * \returns the max length of the peptide_constraint
 */
int get_peptide_constraint_max_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->max_length;
}


/**
 * sets the num_mis_cleavage of the peptide_constraint
 */
void set_peptide_constraint_num_mis_cleavage(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int num_mis_cleavage ///< The maximum mis cleavage of the peptide -in
  )
{
  peptide_constraint->num_mis_cleavage = num_mis_cleavage;
}

/**
 * \returns the num_mis_cleavage of the peptide_constraint
 */
int get_peptide_constraint_num_mis_cleavage(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->num_mis_cleavage;
}

/**
 * sets the mass type of the peptide_constraint
 */
void set_peptide_constraint_mass_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  MASS_TYPE_T mass_type ///< the mass_type for the constraint -in
  )
{
  peptide_constraint->mass_type = mass_type;
}

/**
 * \returns the mass type of the mass_constraint
 */
MASS_TYPE_T get_peptide_constraint_mass_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->mass_type;
}
