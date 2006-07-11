/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.14 $
 * \brief: Object for representing a single peptide.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "objects.h"
#include "mass.h"
#include "peptide.h"
#include "protein.h"


/**
 * \struct peptide
 * \brief A subsequence of a protein.
 */
struct peptide {
  char* sequence;       ///< A pointer to the peptide sequence.
  unsigned char length; ///< The length of the peptide
  float peptide_mass;   ///< The peptide's mass.
};

/**
 * \struct peptide_constraint
 * \brief Object to represent constraints which a peptide may or may not
 *  satisfy.
 * def TRYPTIC: a protein that ends with either K or R and 
 *              any other K and R in the sequence must be followed by a P
 */
struct peptide_constraint {
  BOOLEAN_T is_tryptic; ///< The type of peptides, is there a possibility to be TRYPTIC.
  float min_mass; ///< The minimum mass of the peptide
  float max_mass; ///< The maximum mass of the peptide
  int min_length; ///< The minimum length of the peptide
  int max_length; ///< The maximum length of the peptide
};
 
/**
 * \struct residue_iterator
 * \brief Object to iterate over the residues in a peptide, starting at the
 * first residue of the peptide, and proceeding in order.
 */
struct residue_iterator {
  PEPTIDE_T*  peptide; ///< The peptide whose residues to iterate over.
  int     residue_idx; ///< The index of the current peak
};

/**
 * \returns An (empty) peptide object.
 */
PEPTIDE_T* allocate_peptide(void){
  (PEPTIDE_T*) peptide = (PEPTIDE_T*)mycalloc(1, sizeof(PEPTIDE_T));
  return peptide;
}

/**
 * \returns The mass of the given peptide.
 */
float calc_peptide_mass(
  PEPTIDE_T* peptide ///< the query peptide -in
  )
{
  float peptide_mass = 0;
  RESIDUE_ITERATOR_T * residue_iterator = new_residue_iterator(peptide);
  
  while(residue_iterator_has_next(residue_iterator)){
    peptide_mass += get_amino_acid_mass(residue_iterator_next(residue_iterator));
  }
  return peptide_mass;
}

/**
 * \returns The mass of the given amino acid.
 */
float get_amino_acid_mass(
  char amino_acid ///< the query amino acid -in
  )
{
 
}

/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
PEPTIDE_T* new_peptide(
  char* my_sequence,        ///< The sequence of the protein that that contains the peptide. -in
  unsigned char length,     ///< The length of the peptide -in
  float peptide_mass       ///< The neutral mass of the peptide -in
  )
{
  PEPTIDE_T* peptide = allocate_peptide();
  set_peptide_sequence( peptide, my_sequence);
  set_peptide_length( peptide, length);
  set_peptide_peptide_mass( peptide, peptide_mass);
  return peptide;
}
  

/** 
 * \returns the neutral mass of the peptide
 */
float get_peptide_neutral_mass(
  PEPTIDE_T* peptide ///< the query peptide -in
  )
{
  return get_peptide_peptide_mass(peptide);
}

/** 
 * \returns the mass of the peptide if it had charge "charge"
 */
float get_peptide_charged_mass(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< charge of peptide -in
    )
{
  return get_peptide_mz(peptide, charge) * charge;
}

//TESTME check if this returns a float
/** 
 * \returns the m/z of the peptide if it had charge "charge"
 */
float get_peptide_mz(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< the charge of peptide -in
    )
{
  return ((get_peptide_peptide_mass(peptide) + MASS_H * charge)/ charge);
}

/**
 * Frees an allocated peptide object.
 */
void free_peptide (
  PEPTIDE_T* peptide ///< peptide to free -in
  )
{
  free(peptide->sequence);
  free(peptide);
}

/**
 * Prints a peptide object to file.
 * mass \t peptide-length \t peptide-sequence \n
 */
void print_peptide(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  FILE* file  ///< the out put stream -out
  )
{
  fprintf(file,"%.2f\t",peptide->peptide_mass);
  fprintf(file,"%d\t",peptide->length);
  fprintf(file,"%s\n",peptide->sequence);
}

/**
 * Copies peptide object src to dest.
 * dest must be a heap allocated peptide
 */
void copy_peptide(
  PEPTIDE_T* src, ///< source peptide -in
  PEPTIDE_T* dest ///< destination peptide -out
  )
{
  set_peptide_sequence(dest, get_peptide_sequence(src));
  set_peptide_length(dest, get_peptide_length(src));
  set_peptide_peptide_mass(dest, get_peptide_peptide_mass(src));
}

/**
 * Parses a peptide from file.
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_peptide_file(
  PEPTIDE_T* peptide,
  FILE* file
  )
{
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  char* sequence;
  float peptide_mass;
  int peptide_length;

  line_length =  getline(&new_line, &buf_length, file);
  if( (sscanf(new_line,"%f\t%d\t%s\n", // parse peptide line
              &peptide_mass, &peptide_length, sequence) != 3) && 
      (!isalpha((int)sequence[0]))){
    free(new_line);
    return FALSE;
  }
  set_peptide_sequence(sequence);  // CHECK IF THIS READS the /n as well !!!!!
  set_peptide_peptide_mass(peptide_mass);
  set_peptide_length(peptide_length);
  free(new_line);
  return TRUE;
}

/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* allocate_peptide_constraint(void){
  (PEPTIDE_CONSTRAINT_T*) peptide_constraint =
    (PEPTIDE_CONSTRAINT_T*)mycalloc(1, sizeof(PEPTIDE_CONSTRAINT_T));
  return peptide_constraint;
}

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PEPTIDE_CONSTRAINT_T* new_peptide_constraint(
  BOOLEAN_T is_tryptic, ///< The type of peptides, is there a possibility to be TRYPTIC -in
  float min_mass, ///< the minimum mass -in
  float max_mass, ///< the maximum mass -in
  int min_length, ///< the minimum length of peptide -in
  int max_length  ///< the maximum lenth of peptide -in
  )
{
  (PEPTIDE_CONSTRAINT_T*) peptide_constraint =
    allocate_peptide_constraint();

  set_peptide_constraint_is_tryptic(peptide_constraint, is_tryptic);
  set_peptide_constraint_min_mass(peptide_constraint, min_mass);
  set_peptide_constraint_max_mass(peptide_constraint, max_mass);
  set_peptide_constraint_min_length(peptide_constraint, min_length);
  set_peptide_constraint_max_length(peptide_constraint, max_length);
  return peptide_constraint;
}


//FIXME add is_tryptic?
/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T peptide_constraint_is_satisfied(
    PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraints to enforce -in
    PEPTIDE_T* peptide ///< the query peptide -in
    )
{
  if(get_peptide_length(peptide) <= get_peptide_constraints_max_length(peptide_constraint) &&
     get_peptide_length(peptide) >= get_peptide_constraints_min_length(peptide_constraint) &&
     get_peptide_peptide_mass(peptide) <= get_peptide_constraints_max_mass(peptide_constraint) &&
     get_peptide_peptide_mass(peptide) >= get_peptide_constraints_min_mass(peptide_constraint)
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
  free(peptide_constraint);
}

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
)
{
  free(peptide->sequence);
  int sequence_length = strlen(sequence) +1; //+\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);

  peptide->sequence =
    strncpy(copy_sequence, sequence, sequence_length);  
}

/**
 * \returns the sequence of peptide
 * returns a char* to a heap allocated copy of the sequence
 * user must free the memory
 */
char* get_peptide_sequence(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 )
{
  int sequence_length = strlen(peptide->sequence) +1; //+\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);
  return strncpy(copy_sequence, peptide->sequence, sequence_length);  
}

/**
 * sets the sequence length of the peptide
 */
void set_peptide_length( 
  PEPTIDE_T* peptide,  ///< the peptide to set the length -out
  unsigned char length  ///< the length of sequence -in
  )
{
  peptide->length = length;
}

/**
 *\returns the sequence length of the peptide
 */
unsigned char get_peptide_length( 
  PEPTIDE_T* peptide  ///< the peptide to query the length -in
  )
{
  return peptide->length;
}

/**
 * sets the peptide mass
 */
void set_peptide_peptide_mass( 
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  float peptide_mass  ///< the mass of the peptide - in
  )
{
  peptide->peptide_mass = peptide_mass;
}

/**
 * \returns the peptide mass
 */
float get_peptide_peptide_mass( 
  PEPTIDE_T* peptide  ///< the peptide to query the mass -in
  )
{
  return peptide->peptide_mass;
}



/**
 * sets the peptide type of the peptide_constraint
 */
void set_peptide_constraint_is_tryptic(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out
  BOOLEAN_T is_tryptic ///< the type of the peptide constraint - in
  )
{
  peptide_constraint->is_tryptic = is_tryptic;
}

/**
 * \returns the peptide type of the peptide_constraint
 */
BOOLEAN_T get_peptide_constraint_is_tryptic(
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraint to query -in
  )
{
  return peptide_constraint->is_tryptic;
}

/**
 * sets the min mass of the peptide_constraint
 */
void set_peptide_constraint_min_mass(
  PEPTIDE_CONSTRAINT_T* peptide_constraint, ///< the peptide constraint to set -out
  float min_mass  ///< the min mass of the peptide constraint - in
  )
{
  peptide_constraint->min_mass = min_mass;
}

/**
 * \returns the min mass of the peptide_constraint
 */
float get_peptide_constraint_min_mass(
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
  float max_mass  ///< the max mass of the peptide constraint - in
  )
{
  peptide_constraint->max_mass = max_mass;
}

/**
 * \returns the max mass of the peptide_constraint
 */
float get_peptide_constraint_max_mass(
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
 */
void set_peptide_constraint_max_length(
  PEPTIDE_CONSTRAINT_T* peptide_constraint,///< the peptide constraint to set -out 
  int max_length  ///< the max length of the peptide constraint - in
  )
{
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
 * Iterator
 */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  PEPTIDE_T* peptide ///< peptide sequence to iterate -in
  )
{
  (RESIDUE_ITERATOR_T*) residue_iterator =
    (RESIDUE_ITERATOR_T*)mycalloc(1, sizeof(RESIDUE_ITERATOR_T));
  
  residue_iterator->peptide =  peptide;
  residue_iterator->residue_idx = 0;
  return residue_iterator;
}        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  )
{
  free(residue_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional residues to iterate over, FALSE if not.
 */
BOOLEAN_T residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  )
{
  return (residue_iterator->residue_idx < residue_iterator->peptide->length);
}

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  )
{
  ++residue_idx;
  return residue_iterator->peptide->sequence[residue_idx-1];
}


  
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

