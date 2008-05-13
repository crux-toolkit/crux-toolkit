/**
 * \file modifications.c
 * \brief Datatypes and methods for amino acid modifications
 *
 * Two data structures define modifications.  The AA_MOD_T is the most
 * basic type.  It is the information provided by the user: mass
 * change caused by this mod, amino acids which may be modified in
 * this way, and the maximum number of this type of modification which
 * may occur on one peptide.  AA_MODs are defined here.
 * A collection of AA_MODs that may occur
 * on some peptide are represented as a PEPTIDE_MOD_T.  This stores
 * a list of AA_MODS and the net mass change experienced by the
 * peptide.  PEPTIDE_MODs are defined in peptide_modifications.h
 * AA_MODs are instantiated once after parsing the parameter file.  All
 * possible PEPTIDE_MODs are calcualted once and reused for each
 * spectrum search.  One PEPTIDE_MOD corresponds to one mass window
 * that must be searched.
 * 
 * $Revision: 1.1.2.8 $
 */

#include "modifications.h"

/* Private constants */
//enum { MAX_PROTEIN_SEQ_LENGTH = 40000 };

/* Private global variables */
char mod_sqt_symbols[MAX_AA_MODS] = {'*', '@', '#', '^', '~', '%', 
                                     '$', '&', '!', '?', '+'};
// FIXME: these need to be changed to bitmasks
int mod_id_masks[MAX_AA_MODS] = {1,2,3,4,5,6,7,8,9,10,11};
//MOD_SEQ_NULL = (MODIFIED_AA_T)('Z' - 'A' + 1); 

/* Private data types, typedefed in objects.h */

/**
 * \struct _aa_mod
 * 
 *  Modification at the amino acid level.  A single mass change that can
 *  occur on any of the residues listed.  This information is given by
 *  the user in the parameter file.  Also stores a character symbol
 *  assigned at runtime to be used in the sqt result file and an
 *  integer bitmask to be used to give each aa_mod a unique identifier.
 */
struct _aa_mod{ 
  double mass_change;  ///< the amount by which the mass of the residue changes
  BOOLEAN_T* aa_list;  ///< an array indexed by AA, true if can be modified
  int max_per_peptide; ///< the maximum number of mods per peptide
  MOD_POSITION_T position; ///< where the mod can occur in the pep/prot
  int max_distance;        ///< the max distance from the protein terminus
  char symbol;         ///< the character to represent the mod in sqt files
  int identifier;      ///< the offset/bitmask assigned to this mod for unique
                       //identification, used with MODIFIED_AA
};

/* Definitions of public methods */

/**
 * \brief Allocate an AA_MOD, including the aa_list and initialize all
 * fields to default values.  Symbol and unique identifier are set
 * according to index.  
 * \returns A heap allocated AA_MOD with default values.
 */
AA_MOD_T* new_aa_mod(int mod_idx){

  assert(mod_idx > -1 && mod_idx < MAX_AA_MODS );

  AA_MOD_T* mod = (AA_MOD_T*)mymalloc(sizeof(AA_MOD_T));

  mod->mass_change = 0;
  mod->max_per_peptide = 0;
  mod->position = ANY_POSITION;
  mod->max_distance = MAX_PROTEIN_SEQ_LENGTH;
  mod->symbol = mod_sqt_symbols[mod_idx];
  mod->identifier = mod_id_masks[mod_idx];
  
  // allocate the aa lists for mods 
  mod->aa_list =                        // all 0's?
    (BOOLEAN_T*)mycalloc(AA_LIST_LENGTH, sizeof(AA_MOD_T)); 
  
  // initialize to FALSE
  int aa_idx = 0;
  for(aa_idx=0; aa_idx < AA_LIST_LENGTH; aa_idx++){
    mod->aa_list[aa_idx] = FALSE;
  }
  return mod;
}

/**
 * \brief Free the memory for an AA_MOD including the aa_list.
 */
void free_aa_mod(AA_MOD_T* mod){
  if( mod ){
    if( mod->aa_list ){free(mod->aa_list);}
    free(mod);
  }
}

//FIXME: implment these
/**
 * \brief Converts a MODIFIED_AA into a char, effectively unmodifying it.
 * \returns The unmodified char representation of an aa.
 */
char modified_aa_to_char(MODIFIED_AA_T aa){
  return (char)aa;
}

/**
 * \brief Converts a char representation of an aa to an unmodified aa
 * of type MODIFIED_AA_T.  Requires 'A' <= char <= 'Z'. 
 * \returns The MODIFIED_AA_T represnetation of an aa.
 */
MODIFIED_AA_T char_aa_to_modified(char aa){
  assert( aa >= 'A' && aa <= 'Z' );
  return (MODIFIED_AA_T)(aa - 'A');
}

/**
 * \brief Allocates an array of MODIFIED_AA_T's the same length as
 * sequence and populates it with the MODIFIED_AA_T value that
 * corresponds to each sequence char value.  No modifications are
 * applied to the new array.
 *
 * \returns A newly allocated copy of the sequnce converted to type
 * MODIFIED_AA_T. 
 */
MODIFIED_AA_T* convert_to_mod_aa_seq(char* sequence){

  if( sequence == NULL ){
    carp(CARP_ERROR, "Cannot convert NULL sequence to modifiable characters"); 
    return NULL;
  }

  int seq_len = strlen(sequence);
  MODIFIED_AA_T* new_string = mycalloc( seq_len+1, sizeof(MODIFIED_AA_T) );

  unsigned int seq_idx = 0;
  //  while( sequence[seq_idx] != '\0' ){
  for(seq_idx = 0; seq_idx < strlen(sequence); seq_idx++){
    new_string[seq_idx] = char_aa_to_modified( sequence[seq_idx] );
  }

  // null terminate
  // might need to terminate with something else
  new_string[seq_idx] = MOD_SEQ_NULL;

  return new_string;
}

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq( MODIFIED_AA_T* source){
  if( source == NULL ){
    carp(CARP_ERROR, "Cannot copy NULL sequence of modified_aa's.");
    return NULL;
  }

  // get the length of the source seq
  int i=0;
  while( source[i] != MOD_SEQ_NULL ){
    i++;
  }

  MODIFIED_AA_T* new_seq = mycalloc( i, sizeof(MODIFIED_AA_T) );
  memcpy( new_seq, source, i * sizeof(MODIFIED_AA_T));

  return new_seq;
}

// FIXME: implement this
BOOLEAN_T is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod){
  /*
  int id = mod->identifier;
  if( (aa && id) != 0 ){
    return TRUE;
  } 
  */
  // so it compiles
  carp(110, "%d %d", (int)aa, mod->max_per_peptide);
  // dummy for now
  return FALSE;
 }

/**
 * \brief Determine if this modified amino acid can be modified by
 * this modification.
 *
 * Checks the mod list of modifiable residues to see if aa is in the
 * list.  Also checks to see if aa has already been modified by this
 * mod.  
 * \returns TRUE if it can be modified, else FALSE
 */
BOOLEAN_T is_aa_modifiable
 (MODIFIED_AA_T aa, ///< The sequence amino acid to be modified
  AA_MOD_T* mod){   ///< The aa_mod to modify it with

  if( is_aa_modified(aa, mod) == TRUE ){
    return FALSE;
  }
  //  if( mod->aa_list[ modified_aa_to_char(aa) - 'A' ] == TRUE ){
  if( mod->aa_list[ (int)modified_aa_to_char(aa) ] == TRUE 
      && ! is_aa_modified(aa, mod)){
    return TRUE;
  }
  // else not in list or already modified by this mod
  return FALSE;
}

/**
 * \brief Adds a modification to a MODIFIED_AA_T.
 *
 * Assumes that the aa is modifiable, no explicit check.  If the aa is
 * already modified for the mod, no change to aa.
 */
void modify_aa(MODIFIED_AA_T* aa, AA_MOD_T* mod){
  if( aa == NULL || mod == NULL ){
    carp(CARP_ERROR, "Cannot modify aa.  Either aa or mod NULL.");
    return;
  }
}

/**
 * \brief Finds the list of modifications made to an amino acid.
 *
 * Allocates a list of length(possible_mods) and fills it with pointers
 * to the modifications made to this aa as defined by
 * is_aa_modified().  Returns 0 and sets mod_list to NULL if the amino
 * acid is unmodified
 *
 * \returns The number of modifications made to this amino acid.
 */
/*
int get_aa_mods(MODIFIED_AA_T aa, 
                int num_possible,
                AA_MOD_T* possible_mods, 
                AA_MOD_T*** mod_list){
*/
  /*
    int count = 0;
    mod_list* = new AA_MOD_T[num_possible]
    for 0 to num_possible
      if( is_aa_modified() ){
        mod_list[count] = this mod
        count += 1
      endif
    last num_possible
    if( count == 0 )
       free mod_list
    return count;
   */
//}


/********
 Existing methods to change

 FILE match_collection.c:

 print_sqt_header()
 serialize_headers()
 print_match_collection_sqt()  Actually, the change may only occur in peptide
 serialize_psm_features()
 extend_match_collection()    This is what reads psm header file info

 FILE peptide.c:

 serialize_peptide() minimum: number of mods, two lists--one of
                              identifiers and one of sequnce indexes
 parse_peptide()     if modified, fill sequence, add id to indexed aa


 FILE mass.c:
 get_mass_amino_acid()  change char to short?, new method?

QUESTIONS
* do we require that the parameter file with the mods be included for
  crux_analyae_matches?

* is there a slicker way to index char aa's other than aa - 'A' ?

create branch, add these files to branch
find cvs branching tutorial, send to AK
 */

/**
 * print all fields in aa mod.  For debugging
 */
void print_a_mod(AA_MOD_T* mod){
  printf("AMOD: mass %.2f, max per %d, max dist %d, symb %c, aa list ",
         mod->mass_change, mod->max_per_peptide, mod->max_distance,
         mod->symbol);

  int i = 0;
  for(i=0; i<AA_LIST_LENGTH; i++){
    if( mod->aa_list[i] == TRUE ){
      printf("%c", (i + 'A'));
    }
  }
  printf("\n");
}

/* Setters and Getters */

/**
 * \brief Set the mass change caused by this modification.
 * \returns void
 */
void aa_mod_set_mass_change(AA_MOD_T* mod, double mass_change){
  mod->mass_change = mass_change;
}

/**
 * \brief Get the mass change caused by this modification.
 * \returns The mass change caused by this modification.
 */
double aa_mod_get_mass_change(AA_MOD_T* mod){
  return mod->mass_change;
}

/**
 * \brief Access to the aa_list of the AA_MOD_T struct.  This pointer
 * can be used to get or set the list of residues on which this mod
 * can be placed.
 * \returns A pointer to the list of amino acids on which this mod can
 * be placed.
 */
BOOLEAN_T* aa_mod_get_aa_list(AA_MOD_T* mod){
  return mod->aa_list;
}

/**
 * \brief Set the maximum number of times this modification can be
 * placed on one peptide.
 * \returns void
 */
void aa_mod_set_max_per_peptide(AA_MOD_T* mod, int max){
  mod->max_per_peptide = max;
}

/**
 * \brief Get the maximum number of times this modification can be
 * placed on one peptide.  
 * \returns The max times per peptide this mod can be placed.
 */
int aa_mod_get_max_per_peptide(AA_MOD_T* mod){
  return mod->max_per_peptide;
}

/**
 * \brief Set the maximum distance from the protein terminus that the
 * modification can be placed.  Which terminus (C or N) is determined
 * by the position type.  To indicate no position restriction, set to
 * MAX_PROTEIN_SEQ_LENGTH. 
 * \returns void
 */
void aa_mod_set_max_distance(AA_MOD_T* mod, int distance){
  mod->max_distance = (distance == -1) ? MAX_PROTEIN_SEQ_LENGTH : 
                                         distance;
}

/**
 * \brief Get the maximum distance from the protein end that the
 * modification can be placed.  Will be MAX_PROTEIN_SEQ_LENGTH if
 * position type is ANY_POSITION.
 * \returns Maximum distance from protein terminus at which mod can be
 * placed. 
 */
int aa_mod_get_max_distance(AA_MOD_T* mod){
  return mod->max_distance;
}

/**
 * \brief Set the position type of an aa_mod.
 * \returns void
 */
void aa_mod_set_position(AA_MOD_T* mod, MOD_POSITION_T position){
  mod->position = position;
}

/**
 * \brief Where in the peptide can the modification be placed.
 * \returns ANY_POSITION for standard mods; C_TERM or N_TERM for those
 * that can only be placed at the ends of the peptide.
 */
MOD_POSITION_T aa_mod_get_position(AA_MOD_T* mod){
  return mod->position;
}

/**
 * \brief The character used to uniquely identify the mod in the sqt file.
 * \returns The character identifier.
 */
char aa_mod_get_symbol(AA_MOD_T* mod){
  return mod->symbol;
}

/**
 * \brief The bitmask used to uniquely identify the mod.
 * \returns The short int bitmask used to identify the mod.
 */
int aa_mod_get_identifier(AA_MOD_T* mod){
  return mod->identifier;
}


