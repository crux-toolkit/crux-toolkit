/**
 * \file peptide_src.h
 * $Revision: 1.12 $
 * \brief Object for mapping a peptide to it's parent protein.
 */
#ifndef PEPTIDE_SRC_H
#define PEPTIDE_SRC_H
#include <stdio.h>

#include "utils.h"
#include "mass.h"

#include "objects.h"
#include "carp.h"
#include "peptide_constraint.h"

/**
 * \returns An (empty) peptide_src object.
 */
PEPTIDE_SRC_T* allocate_peptide_src(void);

/**
 * \returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user
 * specified parameters 
 */
PEPTIDE_SRC_T* new_peptide_src(
                               //  PEPTIDE_TYPE_T peptide_type, 
    ///< the peptide type for the corresponding protein -in
DIGEST_T digest,
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< peptide start index in protein sequence, first is 1 -in
  );

/**
 *\returns an array of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for  peptide is known
 */
PEPTIDE_SRC_T* new_peptide_src_array(
  int size ///< the size of the peptide_src array -in
  );

/**
 * \brief Fill in the values from the original array into the new array.
 * Assumes that the new array has been allocated by new_peptide_src_array().
 */
void copy_peptide_src_array(PEPTIDE_SRC_T* original_array, 
                            PEPTIDE_SRC_T* new_array, 
                            int array_size);

/**
 *\returns a linklist of PROTEIN_PEPTIDE_SRC object
 * only used in index.c, when the peptide src count for peptide is known
 */
PEPTIDE_SRC_T* new_peptide_src_linklist(
  int size ///< the size of the peptide_src array -in
  );

/**
 *\returns the PROTEIN_PEPTIDE_SRC object in the array with the index
 * index starts at 0.
 * only used in index.c, when the peptide src count for  peptide is known
 */
void set_peptide_src_array(
  PEPTIDE_SRC_T* src_array , ///< the working peptide src_arry -out
  int array_idx, ///< array index of the peptide_src to set
  //PEPTIDE_TYPE_T peptide_type, ///< the peptide type for the corresponding protein -in
DIGEST_T digest,
  PROTEIN_T* parent_protein, ///< the parent of this preptide -in
  int start_idx ///< start index of the peptide in the protein sequence -in
  );

/**
 * Frees the entire allocated peptide_src linklist object
 * Assumes that peptide src is Link list implementation
 */
void free_peptide_src(
  PEPTIDE_SRC_T* peptide_src ///< object to free -in 
  );

/**
 * Frees the an individual allocated peptide_src object
 * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
 */
void free_one_peptide_src(
  PEPTIDE_SRC_T* peptide_src///< object to free -in 
  );

/**
 * Prints a peptide object to file.
 */
void print_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< object to print -in 
  FILE* file  ///< the out put stream -out
  );

/**
 * \brief Read in the peptide_src objects from the given file and
 * assosiated them with the given peptide.  
 * Proteins for the pepitde_src are found in the given database.  If
 * database is NULL, does not set proteins.  (This option is used for
 * sorting index files while creating index.)  Either array or 
 * linked list implementation of multiple peptide_src is used based on
 * the value of use_array.
 *
 * \returns TRUE if peptide_src's were successfully parsed, else
 * returns FALSE.
 */
BOOLEAN_T parse_peptide_src_tab_delimited(
  PEPTIDE_T* peptide,   ///< assign peptide_src(s) to this peptide
  MatchFileReader& file,           ///< file to read from
  DATABASE_T* database, ///< database containing proteins
  BOOLEAN_T use_array); ///< use array implementation vs. linked list

/**
 * \brief Read in the peptide_src objects from the given file and
 * assosiated them with the given peptide.  
 *
 * Proteins for the pepitde_src are found in the given database.  If
 * database is NULL, does not set proteins.  (This option used for
 * sorting index files while creating index.) Either array or linked
 * list implementation of multiple peptide_src is used based on the
 * value of use_array. 
 *
 * \returns TRUE if peptide_src's were successfully parsed, else
 * returns FALSE and sets peptide's peptide_src member variable to
 * NULL. 
 */
BOOLEAN_T parse_peptide_src(
  PEPTIDE_T* peptide,   ///< assign peptide_src(s) to this peptide
  FILE* file,           ///< file to read from
  DATABASE_T* database, ///< database containing proteins
  BOOLEAN_T use_array);///< use array implementation vs. linked list

/**
 * Copies the entire linklist of peptide_src object src to dest.
 * dest must be a heap allocated peptide_src
 */
void copy_peptide_src(
  PEPTIDE_SRC_T* src, ///< source peptide_src -in
  PEPTIDE_SRC_T* dest ///< destination peptide_src -out
  );

/**
 * sets the peptide type
 * peptide type: TRYPTIC, PARTIALLY_TRYPTIC, NON_TRYPTIC
 */
/*
void set_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PEPTIDE_TYPE_T peptide_type ///< the type of the peptide -in
  );
*/
/**
 * \returns the peptide type with association to the parent protein
 * peptide type: TRYPTIC, N_TRYPTIC, C_TRYPTIC, NON_TRYPTIC
 */
/*
PEPTIDE_TYPE_T get_peptide_src_peptide_type( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );
*/

/**
 * sets the level of digestion
 */
void set_peptide_src_digest( 
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  DIGEST_T digest ///< the type of the peptide -in
  );

/**
 * \returns the level of digestion
 */
DIGEST_T get_peptide_src_digest( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the parent protein
 */
void set_peptide_src_parent_protein(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  PROTEIN_T* parent_protein ///< the parent of this preptide -in  
  );

/**
 * \returns a pointer to the parent protein
 */
PROTEIN_T* get_peptide_src_parent_protein( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the start index of the peptide in the protein sequence
 */
void set_peptide_src_start_idx(
  PEPTIDE_SRC_T* new_association, ///< the peptide_src to set -out   
  int start_idx ///< start index of the peptide in the protein sequence -in
  );

/**
 * \returns the start index of the peptide in the protein sequence
 */
int get_peptide_src_start_idx( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * sets the next peptide_src on the link list
 * assumes that the src_association's next_association feild is NULL
 */
void set_peptide_src_next_association(
  PEPTIDE_SRC_T* src_association, ///< the peptide_src to set -out   
  PEPTIDE_SRC_T* new_association ///< the new peptide_src to add -in   
  );

/**
 * \returns the next peptide_src on the link list
 */
PEPTIDE_SRC_T* get_peptide_src_next_association( 
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 * \returns a pointer to the start of the peptide with in it's parent protein sequence
 */
char* get_peptide_src_sequence_pointer(
  PEPTIDE_SRC_T* peptide_src ///< the query peptide_src -in   
  );

/**
 *\returns the peptide_src strct size, value of sizeof function
 */
int get_peptide_src_sizeof(void);

/**
 * serialize peptide src in binary
 */
void serialize_peptide_src(
  PEPTIDE_SRC_T* peptide_src, ///< peptide_src to serialize -in   
  FILE* file  ///< output file -in   
  );

/**
 * Return the number of bytes taken up by one peptide_src when
 * serialized to file.  Used for skipping past peptide_src in an index
 * file. 
 */
int size_of_serialized_peptide_src();

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
