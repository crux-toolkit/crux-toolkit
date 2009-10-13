/**
 * \file peptide.h 
 * $Revision: 1.52 $
 * \brief Object for representing one peptide.
 */
#ifndef PEPTIDE_H 
#define PEPTIDE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "crux-utils.h"
#include "hash.h"
#include "mass.h"
#include "protein.h"
#include "objects.h"
#include "carp.h"
#include "peptide_constraint.h"
#include "database.h"
#include "modifications.h"
#include "peptide_modifications.h"

//these may be elsewhere
#define MAX_PEPTIDE_LENGTH 255

#ifdef __cplusplus
extern "C" {
#endif


/*  Allocators/deallocators  */

/**
 * \returns An (empty) peptide object.
 */
PEPTIDE_T* allocate_peptide(void);

/**
 *\returns the protein struct size, value of sizeof function
 */
int get_peptide_sizeof(void);

/**
 * \returns A new peptide object, populated with the user specified
 * parameters.
 */
PEPTIDE_T* new_peptide(
  unsigned char length,     ///< The length of the peptide -in
  FLOAT_T peptide_mass,       ///< The neutral mass of the peptide -in
  PROTEIN_T* parent_protein, ///< The parent_protein of this peptide -in
  int start_idx ///< Start index of peptide in the protein sequence -in
  //PEPTIDE_TYPE_T peptide_type ///<  The type of cleavage(TRYPTIC, etc)
  );

/**
 * \brief Allocates a new peptide giving it the values of the source
 * peptide.
 * \returns A newly allocated peptide identical to the source.
 */
PEPTIDE_T* copy_peptide(
  PEPTIDE_T* src ///< source peptide -in
);

/**
 * Merge to identical peptides, copy all peptide_src into one of the peptide
 * peptide_dest, peptide_bye must have at least one peptide src
 * frees the peptide_bye, once the peptide_src are re-linked to the peptide_dest
 * Assumes that both peptides use linklist implemenation for peptide_src
 * \returns TRUE if merge is successful else FALSE
 */
BOOLEAN_T merge_peptides(
  PEPTIDE_T* peptide_dest,
  PEPTIDE_T* peptide_bye
  );
                           
/**
 * Merges two identical peptides by adding the peptide_src of the
 * second to the first.  The second peptide remains unchanged.
 * Does not comfirm identity of peptides.
 * \returns TRUE if merge is successfull.
 */
BOOLEAN_T merge_peptides_copy_src(PEPTIDE_T* peptide_dest,
                                  PEPTIDE_T* peptide_giver);

/**
 * Frees an allocated peptide object.
 * Depending on peptide_src implementation determines how to free srcs
 * This decision is made by global variable PEPTIDE_SRC_USE_LINK_LIST
 */
void free_peptide (
  PEPTIDE_T* peptide ///< peptide to free -in
  );

/*  Getters and Setters  */

/*  Get-set:  mass */
/**
 * sets the peptide src implementation in the peptide object
 * This should be set only once and not be altered
 */
void set_peptide_src_implementation(
  BOOLEAN_T use_link_list ///< does the peptide use link list peptide src
  );

/**
 * sets the peptide mass
 */
void set_peptide_peptide_mass( 
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  FLOAT_T peptide_mass  ///< the mass of the peptide - in
  );

/**
 * \returns the peptide mass
 */
inline FLOAT_T get_peptide_peptide_mass( 
  PEPTIDE_T* peptide  ///< the peptide to query the mass -in
  );

/** 
 * \returns the neutral mass of the peptide
 */
FLOAT_T get_peptide_neutral_mass(
  PEPTIDE_T* peptide ///< the query peptide -in
  );

/** 
 * \returns the mass of the peptide if it had charge "charge"
 */
FLOAT_T get_peptide_charged_mass(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< charge of peptide -in
    );

/** 
 * \returns the m/z of the peptide if it had charge "charge"
 */
FLOAT_T get_peptide_mz(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< the charge of peptide -in
    );

/*  Get-set:  source */

/**
 * sets the peptide_src field in the peptide
 * must pass on a heap allocated peptide_src object
 * does not copy in the object, just the pointer to the object.
 */
void set_peptide_peptide_src(
  PEPTIDE_T* peptide,  ///< the peptide to set -out 
  PEPTIDE_SRC_T* new_association ///< new peptide_src -in
);

/**
 * this method adds the new_association to the end of the existing peptide's 
 * linklist of peptide_srcs
 * must pass on a heap allocated peptide_src object
 * does not copy in the object, just the pointer to the object.
 */
void add_peptide_peptide_src(
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  PEPTIDE_SRC_T* new_association ///< new peptide_src -in
  );

/**
 * this method adds the peptide src array to an EMPTY peptide
 * only used in index.c, when the peptide src count for  peptide is known
 * Any existing peptide_src will lose it's reference
 */
void add_peptide_peptide_src_array(
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  PEPTIDE_SRC_T* peptide_src_array ///< new peptide_src -in
  );

/**
 * returns a point to the peptide_protein_association field of the peptide
 */
PEPTIDE_SRC_T* get_peptide_peptide_src(
  PEPTIDE_T* peptide  ///< the peptide to query the peptide_peptide_src -in
);

/**
 * get the peptide->first peptide_src->parent protein->database
 */
DATABASE_T* get_peptide_first_src_database(
  PEPTIDE_T* peptide ///< working peptide -in
);

/**
 * returns a pointer to the peptide's first parent protein field of the peptide
 */
PROTEIN_T* get_peptide_parent_protein(
  PEPTIDE_T* peptide  ///< the peptide to query the parent_protein -in
  );

/*  Get-set:  sequence */

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
 * \brief Get the sequence of a peptide.
 * Goes to the first peptide_src to gain sequence, thus must have at
 * least one peptide src 
 * \returns A newly allocated copy of the sequence.
 */
char* get_peptide_sequence(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 );

/**
 * \returns a pointer to the start of peptide sequence with in it's protein parent sequence, thus does not have terminating signe until end of parent protein
 * goes to the first peptide_src to find the location of start, thus must have at least one peptide src
 * should not print, will result in printing the entire protein sequence
 */
char* get_peptide_sequence_pointer(
  PEPTIDE_T* peptide ///< peptide to query sequence -in
  );

/**
 * \brief Formats the sequence of the peptide with each flanking AA.
 * 
 * Format is "X.peptide_sequence.X", where "X" is a flanking amino acid.
 * "X", is printed as "-" if there is no flanking sequence.
 * Goes to the first peptide_src to gain sequence, thus must have at
 * least one peptide src 
 * \returns A newly allocated string with the sqt-formated peptide sequence.
 */
char* get_peptide_sequence_sqt(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 );

/**
 * \brief Formats the sequence of the peptide from a particular
 * peptide_src.
 *
 * Is called by get_peptide_sequence_sqt()
 * Format is "X.peptide_sequence.X", where "X" is a flanking amino acid.
 * "X", is printed as "-" if there is no flanking sequence.
 * Goes to the first peptide_src to gain sequence, thus must have at
 * least one peptide src 
 *
 * \returns A newly allocated string with the sqt-formated peptide sequence.
 */
char* get_peptide_sequence_from_peptide_src_sqt(
 PEPTIDE_T* peptide, ///< peptide to query sequence -in
 PEPTIDE_SRC_T* peptide_src ///< peptide_src -in 
 );

/**
 * \brief Return a char for the amino acid c-terminal to the peptide
 * in the peptide src at the given index.
 *
 * \returns A char (A-Z) or - if peptide is the first in the protein.
 */
char get_peptide_c_term_flanking_aa(
 PEPTIDE_T* peptide   ///< peptide of interest
 );

/**
 * \brief Return a char for the amino acid n-terminal to the peptide
 * in the peptide src at the given index.
 *
 * \returns A char (A-Z) or - if peptide is the last in the protein.
 */
char get_peptide_n_term_flanking_aa(
 PEPTIDE_T* peptide   ///< peptide of interest
 );

/**
 * \brief Add a modification to a peptide.
 *
 * Adds the modified sequence to the peptide and changes the peptide
 * mass based on the mass change in the peptide_mod.
 * \returns void
 */
void set_peptide_mod(PEPTIDE_T* peptide,     ///< peptide to be modified
                     MODIFIED_AA_T* mod_seq, ///< modified seq to add
                     PEPTIDE_MOD_T* pep_mod  ///< mod that made the seq
);

/**
 * \brief Get the modified aa sequence
 *
 * If the peptide has no modifications, create a sequence of
 * MODIFIED_AA_T's in which none of them are actually modified.
 * \returns A newly allocated copy of the sequence of MODIFIED_AA_Ts.
 */
//MODIFIED_AA_T* get_peptide_modified_sequence( // why is this not working??!!
unsigned short* get_peptide_modified_aa_sequence(
 PEPTIDE_T* peptide
 );

/**
 * \brief Get the modified aa sequence in string form.
 *
 * If the peptide has no modifications, returns same string as
 * get_peptide_sequence.  If modified, adds the mod symbols to the string.
 * \returns A newly allocated string of the peptide sequence including
 * any modifications.
 */
char* get_peptide_modified_sequence(
 PEPTIDE_T* peptide
 );

/*  Getters requiring calculation */
int count_peptide_modified_aas(PEPTIDE_T* peptide);

/**
 * \returns The mass of the given peptide.
 */
FLOAT_T calc_sequence_mass(
  const char* peptide, ///< the query peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  );

/**
 * \returns The mass of the given peptide.
 */
FLOAT_T calc_peptide_mass(
  PEPTIDE_T* peptide, ///< the query peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  );

/**
 * \returns The hydrophobicity of the given peptide, as in Krokhin (2004).
 */
FLOAT_T calc_krokhin_hydrophobicity(
  PEPTIDE_T* peptide ///< the query peptide -in
  );

/**
 * Examines the peptide sequence and counts how many tryptic missed
 * cleavage sites exist. 
 *\returns the number of missed cleavage sites in the peptide
 */
int get_peptide_missed_cleavage_sites(
  PEPTIDE_T* peptide  ///< the peptide to query -in
);

/**
 * \brief Find the distance from the n-terminus of the source protein
 * to the n-terminus of the peptide.  
 * In the case of multiple source proteins, return the smallest
 * distance.
 * \returns The distance from the protein n-terminus.
 */
int get_peptide_n_distance(PEPTIDE_T* peptide);

/**
 * \brief Find the distance from the c-terminus of the source protein
 * to the c-terminus of the peptide.
 * In the case of multiple source proteins, return the smallest
 * distance.
 * \returns The distance from the protein c-terminus.
 */
int get_peptide_c_distance(PEPTIDE_T* peptide);

/**
 * Creates a heap allocated hash_value for the peptide that should
 * uniquely identify the peptide
 *\returns the string of "<first src protein idx><start idx><length>"
 */
char* get_peptide_hash_value( 
  PEPTIDE_T*  peptide ///< The peptide whose residues to iterate over.
  );

/**
 * \brief Return a randomly shuffled version of the given peptide's 
 * sequence as an array of char (A-Z).  Based on the peptide type,
 * will leave the end(s) unchanged to preserve the tryptic property. 
 * 
 *\returns A newly-allcoated char array of the shuffled sequence.
 */
char* generate_shuffled_sequence(
  PEPTIDE_T* peptide ///< The peptide to shuffle -in 
  //PEPTIDE_TYPE_T peptide_type 
    ///< tryptic status to enforce on the shuffled sequence
  );

/**
 * \brief Return a reversed version of the given peptide's sequence as
 * an array of char (A-Z).  Leave the first and last residue
 * unchanged.  If the reversed sequence is identical to the target,
 * shuffle the sequence instead.
 *
 * \returns A newly-allocated char array of the reversed sequence.
 */
char* generate_reversed_sequence(
  PEPTIDE_T* peptide ///< The peptide to shuffle -in 
  );

/**
 * \brief Return a randomly shuffled version of the given peptide's 
 * sequence as an array of MODIIFIED_AA_T.  Based on the peptide type,
 * will leave the end(s) unchanged to preserve the tryptic property.
 * 
 *\returns A newly-allcoated MODIFIED_AA_T array of the shuffled sequence.
 */
MODIFIED_AA_T* generate_shuffled_mod_sequence(
  PEPTIDE_T* peptide  ///< The peptide to shuffle -in
  //PEPTIDE_TYPE_T peptide_type 
  ///< tryptic status to enforce on the shuffled sequence
  // not currently used
  );

/**
 * \brief Return a reversed version of the given peptide's sequence as
 * an array of MODIFIED_AA_T.  Leave the first and last residue
 * unchanged.  If the reversed sequence is identical to the target,
 * shuffle the sequence instead.
 *
 * \returns A newly-allocated MODIFIED_AA_T array of the reversed sequence.
 */
MODIFIED_AA_T* generate_reversed_mod_sequence(
  PEPTIDE_T* peptide ///< The peptide to shuffle -in 
  );


/*  Comparisons for sorting  */

/**
 * Compare peptide sequence
 * \returns TRUE if peptide sequence is identical else FALSE
 */
BOOLEAN_T compare_peptide_sequence(
  PEPTIDE_T* peptide_one,
  PEPTIDE_T* peptide_two
  );

/**
 * compares two peptides with the lexical sort type
 * for qsort
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_lexical_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  );


/**
 * compares two peptides with the mass sort type
 * if peptide mass is identical sort by lexicographical order
 * used for qsort function
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_mass_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  );


/**
 * compares two peptides with the length sort type
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_length_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  );

/**
 * Compare peptide mass
 * \returns 0 if peptide mass is identical else 1 if peptide_one is larger, -1 if peptide_two is larger
 */
int compare_peptide_mass(
  PEPTIDE_T* peptide_one,
  PEPTIDE_T* peptide_two
  );

/*  Printing / parsing       */

/**
 * Prints a peptide object to file.
 * prints all peptide_src object it's associated 
 * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
 *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
 * prints in correct format for generate_peptide
 */
void print_peptide_in_format(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  BOOLEAN_T flag_out, ///< print peptide sequence? -in
  //BOOLEAN_T trypticity_opt, ///< print trypticity of peptide? -in
  FILE* file  ///< the out put stream -out
  );

/**
 * Prints a peptide object to file.
 * ONLY prints peptide_src that match the peptide_src
 * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 * prints in correct format for generate_peptide
 */
void print_filtered_peptide_in_format(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  BOOLEAN_T flag_out, ///< print peptide sequence? -in
  FILE* file  ///< the out put stream -out
  //PEPTIDE_TYPE_T peptide_type ///< the peptide_type of src to print -in
  );

/**
 * Serialize a peptide to a FILE in binary
 * \returns TRUE if serialization is successful, else FALSE
 *
 * The peptide serialization format looks like this:
 *
 *<PEPTIDE_T: peptide struct><int: number of peptide_src>[<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>]+
 * the bracket peptide src information repeats for the number of peptide src listed before the bracket
 * the protein index is the index of the parent protein in the database DATABASE_T
 *
 */
BOOLEAN_T serialize_peptide(
  PEPTIDE_T* peptide,
  FILE* file
  );
 
/**
 * \brief Read in a peptide from a binary file and return it.
 *
 * Assumes the peptide has been written to file using
 * serialize_peptide().  Allocates memory for the peptide and all of
 * its peptide_src's.  Requires a database so that the protein can be
 * set for each peptide_src.  Returns NULL if eof or if file format
 * appears incorrect.
 *
 * \returns A newly allocated peptide or NULL
 */
PEPTIDE_T* parse_peptide(
  FILE* file, ///< the serialized peptide file -in
  DATABASE_T* database,///< the database containing the peptides -in
  BOOLEAN_T use_array  ///< should I use array peptide_src or link list -in  
  );

/**
 * \brief Read in a peptide from a binary file without reading its
 * peptide_src's.
 *
 * This parsing method is for callers that do not want memory
 * allcoated for every peptide in the file.  Caller allocates memory
 * once, parses peptide, checks values, and returns or keeps looking.
 * To get the peptide_src for this peptide, caller uses 
 * fseek(file, peptide_src_file_location, SEEK_SET);
 * parse_peptide_src(peptide, file, database, use_array);
 *
 * Assumes that the peptide has been written to file using
 * serialize_peptide().  
 * \returns TRUE if peptide was successfully parsed or FALSE if it was
 * not. 
 */
BOOLEAN_T parse_peptide_no_src(
  PEPTIDE_T* peptide, ///< memory already allocated for a peptide 
  FILE* file,       ///< file pointing to a serialized peptide
  long int* pepitde_src_file_location);  // use to seek back to peptide_src

/*  Iterators */

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

/**
 * Protein peptide association Iterator
 */

/**
 * Instantiates a new peptide_src_iterator from a peptide.
 * \returns a PEPTIDE_SRC_T object.
 */
PEPTIDE_SRC_ITERATOR_T* new_peptide_src_iterator(
  PEPTIDE_T* peptide ///< peptide's fields to iterate -in
  );

/**
 * Frees an allocated peptide_src_iterator object.
 */
void free_peptide_src_iterator(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator ///< free this object -in
  );

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptide_srcs to iterate over, FALSE if not.
 */
BOOLEAN_T peptide_src_iterator_has_next(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator///< the query iterator -in
  );

/**
 * \returns The next peptide_srcs in the peptide.
 */
PEPTIDE_SRC_T* peptide_src_iterator_next(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator///< the query iterator -in
  );

/**
 * \brief Builds a comma delimited string listing the protein ids
 * for the sources of a peptide.
 *
 * \returns a pointer to the string. Caller is responsible for freeing memeory.
 * If peptide has no sources returns NULL.
 */
char *get_protein_ids(PEPTIDE_T *peptide);

/**
 * \brief Builds a comma delimited string listing the flanking amino acids
 * for the sources of a peptide.
 *
 * \returns a pointer to the string. Caller is responsible for freeing memeory.
 * If peptide has no sources returns NULL.
 */
char *get_flanking_aas(PEPTIDE_T *peptide);

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
