/**
 * \file peptide.h 
 * $Revision: 1.10 $
 * \brief Object for representing one peptide.
 *****************************************************************************/
#ifndef PEPTIDE_H 
#define PEPTIDE_H

#include "utils.h"
#include "mass.h"
#include <stdio.h>


/**
 * \typedef PEPTIDE_T
 * \brief A peptide subsequence of a protein
 */
typedef struct peptide PEPTIDE_T;

/**
 * \typedef PEPTIDE_CONSTRAINT_T
 * \brief An object representing constraints which a peptide may or may not
 * satisfy.
 */
typedef struct peptide_constraint PEPTIDE_CONSTRAINT_T;

/**
 * \typedef RESIDUE_ITERATOR_T 
 * \brief An object to iterate over the residues in a peptide
 */
typedef struct residue_iterator RESIDUE_ITERATOR_T;

/**
 * \brief The enum for peptide type, with regard to trypticity.
 */
enum _peptide_type { TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC}; 

/**
 * \typedef PEPTIDE_TYPE_T 
 * \brief The typedef for peptide type, with regard to trypticity.
 */
typedef enum _peptide_type PEPTIDE_TYPE_T;


/**
 * \returns The mass of the given peptide.
 */
float calc_peptide_mass(PEPTIDE_T* peptide);


/**
 * \returns An (empty) peptide object.
 */
PEPTIDE_T* allocate_peptide(void);

/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
PEPTIDE_T* new_peptide(
  char* my_sequence,        ///< The sequence of the protein that that contains the peptide.
  unsigned char length,     ///< The length of the peptide
  double peptide_mass       ///< The neutral mass of the peptide
);
  

/** 
 * \returns the neutral mass of the peptide
 */
float get_peptide_neutral_mass(PEPTIDE_T* peptide);

/** 
 * \returns the mass of the peptide if it had charge "charge"
 */
float get_peptide_charged_mass(
    PEPTIDE_T* peptide, 
    int charge);

/** 
 * \returns the m/z of the peptide if it had charge "charge"
 */
float get_peptide_mz(
    PEPTIDE_T* peptide, 
    int charge);

/**
 * Frees an allocated peptide object.
 */
void free_peptide (PEPTIDE_T* peptide);

/**
 * Prints a peptide object to file.
 */
void print_peptide(PEPTIDE_T* peptide, FILE* file);

/**
 * Copies peptide object src to dest.
 */
void copy_peptide(
  PEPTIDE_T* src,
  PEPTIDE_T* dest);

/**
 * Parses a peptide from file.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_peptide_file(
  PEPTIDE_T* peptide,
  FILE* file);

/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* allocate_peptide_constraint(void);

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint(PEPTIDE_TYPE_T* peptide_type,
    float min_mass, float max_mass, int min_length, int max_length);

/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T peptide_constraint_is_satisfied(
    PEPTIDE_CONSTRAINT_T* peptide_constraint, 
    PEPTIDE_T* peptide);

/**
 * Frees an allocated residue_iterator object.
 */
void free_peptide_constraint(PEPTIDE_CONSTRAINT_T* peptide_constraint);

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 * Iterator
 */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(PEPTIDE_T* peptide);        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(RESIDUE_ITERATOR_T* residue_iterator);

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional residues to iterate over, FALSE if not.
 */
BOOLEAN_T residue_iterator_has_next(RESIDUE_ITERATOR_T* residue_iterator);

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(RESIDUE_ITERATOR_T* residue_iterator);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
