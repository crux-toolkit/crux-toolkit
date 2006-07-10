/**
 * \file peptide.h 
 * $Revision: 1.11 $
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
float get_peptide_neutral_mass(
  PEPTIDE_T* peptide ///< the query peptide -in
  );

/** 
 * \returns the mass of the peptide if it had charge "charge"
 */
float get_peptide_charged_mass(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< charge of peptide -in
    );

/** 
 * \returns the m/z of the peptide if it had charge "charge"
 */
float get_peptide_mz(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< the charge of peptide -in
    );

/**
 * Frees an allocated peptide object.
 */
void free_peptide (
  PEPTIDE_T* peptide ///< peptide to free -in
  );

/**
 * Prints a peptide object to file.
 */
void print_peptide(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  FILE* file  ///< the out put stream -out
  );

/**
 * Copies peptide object src to dest.
 */
void copy_peptide(
  PEPTIDE_T* src, ///< source peptide -in
  PEPTIDE_T* dest ///< destination peptide -out
  );

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
PEPTIDE_CONSTRAINT_T* new_peptide_constraint(
  PEPTIDE_TYPE_T* peptide_type, ///< the peptide_type -in
  float min_mass, ///< the minimum mass -in
  float max_mass, ///< the maximum mass -in
  int min_length, ///< the minimum length of peptide -in
  int max_length  ///< the maximum lenth of peptide -in
  );

/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T peptide_constraint_is_satisfied(
   PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraints to enforce -in
   PEPTIDE_T* peptide ///< the query peptide -in
   );

/**
 * Frees an allocated peptide_constraint object.
 */
void free_peptide_constraint(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< object to free -in 
  );

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */



/**
 * sets the sequence of the peptide
 * copies in the sequence into a new heap allocated string 
 */
void set_peptide_sequence( 
  PEPTIDE_T* peptide,  ///< the peptide to set the sequence -out
  char* sequence ///< the sequence to copy -in
  );

/**
 * \returns the sequence of peptide
 * returns a char* to a heap allocated copy of the sequence
 * user must free the memory
 */
char* get_peptide_sequence(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 );

/**
 * sets the sequence length of the peptide
 */
void set_peptide_length( 
  PEPTIDE_T* peptide,  ///< the peptide to set the length -out
  unsigned char length  ///< the length of sequence -in
  );

/**
 *\returns the sequence length of the peptide
 */
unsigned char get_peptide_length( 
  PEPTIDE_T* peptide  ///< the peptide to query the length -in
  );

/**
 * sets the peptide mass
 */
void set_peptide_peptide_mass( 
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  float peptide_mass  ///< the mass of the peptide - in
  );

/**
 * \returns the peptide mass
 */
float get_peptide_peptide_mass( 
  PEPTIDE_T* peptide  ///< the peptide to query the mass -in
  );

/**
 * sets the peptide type of the peptide_constraint
 */
void set_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide constraint - in
  );

/**
 * \returns the peptide type of the peptide_constraint
 */
PEPTIDE_TYPE_T* get_peptide_constraint_peptide_type(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the min mass of the peptide_constraint
 */
void set_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide constraint to set -out
  float min_mass  ///< the min mass of the peptide constraint - in
  );

/**
 * \returns the min mass of the peptide_constraint
 */
float get_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  );

/**
 * sets the max mass of the peptide_constraint
 */
void set_peptide_constraint_max_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  float max_mass  ///< the max mass of the peptide constraint - in
  );

/**
 * \returns the max mass of the peptide_constraint
 */
float get_peptide_constraint_max_mass(
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
 * Iterator
 */



/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  PEPTIDE_T* peptide ///< peptide sequence to iterate -in
  );

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional residues to iterate over, FALSE if not.
 */
BOOLEAN_T residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  );

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
