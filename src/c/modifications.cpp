/**
 * \file modifications.cpp
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
 * $Revision: 1.8 $
 */

#include "modifications.h"

/* Private constants */
//enum { MAX_PROTEIN_SEQ_LENGTH = 40000 };

/* Private global variables */
/*
   Modified sequences need both a textual representation as well as a
   way to be encoded.  A peptide sequence is typically an array of
   char.  A modified sequence will now be an array of MODIFIED_AA_T,
   which are shorts.  Each MODIFIED_AA_T is made up of two parts: the
   5 least-significant bits are used to encode the amino acid, and the
   remaining are used to encode the modifications made to the aa.

   Values 0-26 are used to represent amino acids A-Z.  There is a
   function that will take an MODIFIED_AA_T and return the text
   representation of the amino acid.

   The value of the most significant bit is used to represent
   modification by the first AA_MOD_T in the parameter.c list.  The
   second most significant bit represents the second MOD_AA_T in the
   list and so on.  Note that each AA_MOD_T can be applied only once
   to an amino acid.  A MODIFED_AA_T can have between 0 and 11
   modifications applied to it.  There is a function that will take a
   MODIFIED_AA_T* and the index of an AA_MOD_T and return TRUE if the
   MODIFIED_AA_T has had the AA_MOD_T applied to it.

   For a textual representation of a modiciation, there are 11
   non-alpha-numeric characters chosen.  There is a function that will
   take an MODIFIED_AA_T* and return a pointer to an array of
   characters that begins with the letter of the amino acid and is
   followed by all of the symbols corresponding to all of the
   modifications that have been applied to the amino acid.
 */
char mod_sqt_symbols[MAX_AA_MODS] = //{'*', '@', '#', '^', '~', '%', 
  {'*', '#', '@', '^', '~', '%', // like sequest
                                     '$', '&', '!', '?', '+'};

unsigned short mod_id_masks[MAX_AA_MODS] = 
  {
   0x0020,   // 0000 0000 0010 0000 *
   0x0040,   // 0000 0000 0100 0000 #
   0x0080,   // 0000 0000 1000 0000 @
   0x0100,   // 0000 0001 0000 0000 ^
   0x0200,   // 0000 0010 0000 0000 ~
   0x0400,   // 0000 0100 0000 0000 %
   0x0800,   // 0000 1000 0000 0000 $
   0x1000,   // 0001 0000 0000 0000 &
   0x2000,   // 0010 0000 0000 0000 !
   0x4000,   // 0100 0000 0000 0000 ?
   0x8000    // 1000 0000 0000 0000 +
}; 

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
  BOOLEAN_T aa_list[AA_LIST_LENGTH];
                       ///< an array indexed by AA, true if can be modified
  int max_per_peptide; ///< the maximum number of mods per peptide
  MOD_POSITION_T position; ///< where the mod can occur in the pep/prot
  int max_distance;        ///< the max distance from the protein terminus
  char symbol;         ///< the character to represent the mod in sqt files
  MODIFIED_AA_T identifier; ///< the bitmask assigned for unique ID
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
  /*
  mod->aa_list =                        // all 0's?
    (BOOLEAN_T*)mycalloc(AA_LIST_LENGTH, sizeof(BOOLEAN_T)); 
  */
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
    //if( mod->aa_list ){free(mod->aa_list);}
    free(mod);
  }
}

/**
 * \brief Gives the size of the aa_mod struct.  For serialization
 */
int get_aa_mod_sizeof(){
  return sizeof(AA_MOD_T);
}

/**
 * \brief Converts a MODIFIED_AA into a char, effectively unmodifying it.
 * \returns The unmodified char representation of an aa.
 */
char modified_aa_to_char(MODIFIED_AA_T aa){
  aa = aa & 0x001F;         // 0000 0000 0001 1111
  return (char)aa + 'A';
}

/**
 * \brief Converts a char representation of an aa to an unmodified aa
 * of type MODIFIED_AA_T.  Requires 'A' <= char <= 'Z'. 
 * \returns The MODIFIED_AA_T represnetation of an aa.
 */
MODIFIED_AA_T char_aa_to_modified(char aa){
  assert( aa >= 'A' && aa <= 'Z' );
  //  return (MODIFIED_AA_T)(aa - 'A');
  MODIFIED_AA_T mod_aa = (MODIFIED_AA_T)aa - (MODIFIED_AA_T)'A';
  return mod_aa;
}

/**
 * \brief Converts a MODIFIED_AA_T* to it's textual representation,
 * i.e. a letter followed by between 0 and 11 symbols for the
 * modifications made to the amino acid.
 * \returns A newly allocated char* with amino acid and modifciation
 * symbols. 
 */
char* modified_aa_to_string(MODIFIED_AA_T aa){

  int modified_by = 0;
  int mod_idx = 0;
  AA_MOD_T** mod_list = NULL;
  //int total_mods = get_aa_mod_list(&mod_list);
  int total_mods = get_all_aa_mod_list(&mod_list);
  for(mod_idx = 0; mod_idx< total_mods; mod_idx++){
    if( is_aa_modified(aa, mod_list[mod_idx])){
      modified_by++;
    }
  }
  char* return_string = (char*)mymalloc((modified_by+2)*sizeof(char));
                                                  // 2 = aa + null term
  int return_idx = 0;
  return_string[return_idx] = modified_aa_to_char(aa);
  return_idx++;
  for(mod_idx = 0; mod_idx < total_mods; mod_idx++){
    if( is_aa_modified(aa, mod_list[mod_idx])){
      return_string[return_idx] = aa_mod_get_symbol(mod_list[mod_idx]);
      return_idx++;
    }
  }
  return_string[return_idx] = '\0';
  return return_string;
}

/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the symbol for all
 * applied modifications.
 *
 * Assumes that the array is terminated with MOD_SEQ_NULL.
 * \returns A newly allocated array of characters, a text
 * representation of the modified sequence.
 */
char* modified_aa_string_to_string(MODIFIED_AA_T* aa_string, int length){
  if( aa_string == NULL ){
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }

  AA_MOD_T** global_mod_list = NULL;
  int mod_list_length = get_all_aa_mod_list(&global_mod_list);
  int global_idx = 0;

  // count up amino acids and modifications
  int count = 0;
  int mod_str_idx = 0;
  //  while( aa_string[mod_str_idx] != MOD_SEQ_NULL ){
  for(mod_str_idx=0; mod_str_idx< length; mod_str_idx++){
    count++;      // count the aas
    for(global_idx = 0; global_idx < mod_list_length; global_idx++){
      if( is_aa_modified(aa_string[mod_str_idx], 
                         global_mod_list[global_idx]) ){
        count++;  // count all mods
      }
    }
    //    mod_str_idx++;
  }
  //  carp(CARP_DETAILED_DEBUG, "Found %d modified aas out of %i",
  //     count - mod_str_idx, count);
  //int mod_str_len = mod_str_idx;  // keep track of aa_string length

  // create a buffer to hold all aas and mod symbols
  char* return_string = (char*)mymalloc((count+1)*sizeof(char));
  char* return_str_ptr = return_string;

  // for each mod_aa, get the string representation, copy to return str
  //  for(mod_str_idx = 0; mod_str_idx<mod_str_len; mod_str_idx++){
  for(mod_str_idx = 0; mod_str_idx<length; mod_str_idx++){

    char* cur_mod = modified_aa_to_string( aa_string[mod_str_idx] );
    strcpy( return_str_ptr, cur_mod );
    return_str_ptr += strlen(cur_mod);
    free(cur_mod);
  }

  return return_string;
}

/**
 * \brief Takes an array of MODIFIED_AA_T's and returns an array of
 * char's, one for each amino acid in the sequence.  It DOES NOT
 * include any modification symbols.  Use with caution.
 * \returns A newly allocated char* with ONLY amino acids, all
 * modifications are removed.
 */
char* modified_aa_to_unmodified_string(MODIFIED_AA_T* aa_string, int length){

  if( aa_string == NULL ){
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }

  char* new_string = (char*)mycalloc(length+1, sizeof(char));
  int aa_idx = 0;
  for(aa_idx = 0; aa_idx < length; aa_idx++){
    new_string[aa_idx] = modified_aa_to_char(aa_string[aa_idx]);
  }
  new_string[length] = '\0';

  return new_string;
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
MODIFIED_AA_T* convert_to_mod_aa_seq(const char* sequence){

  if( sequence == NULL ){
    carp(CARP_ERROR, "Cannot convert NULL sequence to modifiable characters"); 
    return NULL;
  }

  int seq_len = strlen(sequence);
  MODIFIED_AA_T* new_string = (MODIFIED_AA_T*)mycalloc( seq_len+1, sizeof(MODIFIED_AA_T) );

  unsigned int seq_idx = 0;
  //  while( sequence[seq_idx] != '\0' ){
  for(seq_idx = 0; seq_idx < strlen(sequence); seq_idx++){
    new_string[seq_idx] = char_aa_to_modified( sequence[seq_idx] );
  }

  // null terminate
  new_string[seq_idx] = MOD_SEQ_NULL;

  return new_string;
}

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(MODIFIED_AA_T* source, int length){
  if( source == NULL ){
    carp(CARP_ERROR, "Cannot copy NULL sequence of modified_aa's.");
    return NULL;
  }

  MODIFIED_AA_T* new_seq = (MODIFIED_AA_T*)mycalloc( length + 1, sizeof(MODIFIED_AA_T) );
  memcpy( new_seq, source, length * sizeof(MODIFIED_AA_T));
  new_seq[length] = MOD_SEQ_NULL;

  return new_seq;
}

/**
 * \brief Determine if an array of MODIFIED_AA_T is a palindrome.  
 * Used by reverse_sequence to avoid returning a reversed sequence
 * that is the same as the target.  Ignores the first and last
 * residues. 
 * \returns TRUE if the reversed sequence would be the same as the
 * forward, otherwise FALSE.
 */
BOOLEAN_T modified_aa_seq_is_palindrome(MODIFIED_AA_T* seq, int length){
  if( seq == NULL ){
    return FALSE;
  }

  int left_idx = 1;    // skip first and last residues
  int right_idx = length - 2;

  while( left_idx < right_idx ){
    if( seq[left_idx] != seq[right_idx] ){
      return FALSE;
    }// else, keep checking
    left_idx++;
    right_idx--;
  }

  // if we got to here, they all matched
  return TRUE;
}


/**
 * \brief Write the given aa mod to file in binary format.  Used for
 * serializing csm file headers.
 *
 * \return TRUE if written to file without error, else false.
 */
BOOLEAN_T serialize_aa_mod(AA_MOD_T* a_mod,
                           FILE* file
){
  if( a_mod == NULL || file == NULL ){
    carp(CARP_ERROR, "Cannot serialize aa_mod to file with NULL inputs");
    return FALSE;
  }

  // write all member variables
  fwrite(a_mod, sizeof(AA_MOD_T), 1, file);

  // then write each aa TRUE/FALSE
  /*
  int aa_idx;
  for(aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
    fwrite((a_mod->aa_list)+aa_idx, sizeof(BOOLEAN_T), 1, file);
  }
  */
  return TRUE;
}

/**
 * \brief Read an aa mod from file in binary format as written by
 * serialize_aa_mod().  Overwrites any data in the passed aa_mod.
 * Used for reading in csm files. If FALSE is returned, the passed
 * a_mod will contain unpredictable values.
 *
 * \returns TRUE if aa_mod successfully read, else FALSE.  
 */
BOOLEAN_T parse_aa_mod(AA_MOD_T* a_mod,
                       FILE* file
){
  if( a_mod == NULL || file == NULL ){
    carp(CARP_ERROR, "Cannot parse aa_mod from file with NULL inputs");
    return FALSE;
  }

  // read in struct
  if( fread(a_mod, sizeof(AA_MOD_T), 1, file) != 1 ){
    carp(CARP_ERROR, "Failed to read AA_MOD.");
  };

  // then read in each modifiable aa (TRUE/FALSE)
  /*
  int aa_idx;
  for(aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++){
    fread((a_mod->aa_list)+aa_idx, sizeof(BOOLEAN_T), 1, file);
  }
  */
  return TRUE;
}

/**
 * \brief Frees memory for an array of MODIFIED_AA_Ts.  Assumes is
 * terminated with the MOD_SEQ_NULL value
 */
/*
void free_mod_aa_seq( MODIFIED_AA_T* seq ){

}
*/

/**
 * \brief Check that the list of peptide_modifications from the file of
 * serialized PSMs matches those in the paramter file.
 *
 * If there was no parameter file or if it did not contain any mods,
 * return FALSE.  If the given mod list does not exactly match the
 * mods read from the parameter file (including the order in which
 * they are listed) return FALSE.  If returning false, print a warning
 * with the lines that should be included in the parameter file.
 *
 * \returns TRUE if the given mods are the same as those from the
 * parameter file.
 */
BOOLEAN_T compare_mods(AA_MOD_T** psm_file_mod_list, int file_num_mods){
  AA_MOD_T** mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&mod_list);

  if( num_mods != file_num_mods ){
    return FALSE;
  }
  int mod_idx = 0;
  for(mod_idx=0; mod_idx<num_mods; mod_idx++){
    if( ! compare_two_mods(mod_list[mod_idx], psm_file_mod_list[mod_idx]) ){
      printf("param mod %d is: ", mod_idx);
      print_a_mod(mod_list[mod_idx]);
      printf("file mod %d is: ", mod_idx);
      print_a_mod(psm_file_mod_list[mod_idx]);
      return FALSE;
    }
  }
  return TRUE;
}

/**
 * \brief Compare two mods to see if they are the same, i.e. same mass
 * change, unique identifier, position
 */
BOOLEAN_T compare_two_mods(AA_MOD_T* mod1, AA_MOD_T* mod2){
  if(mod1->mass_change != mod2->mass_change ){
    return FALSE;
  }
  if(mod1->max_per_peptide != mod2->max_per_peptide ){
    return FALSE;
  }
  if(mod1->position != mod2->position ){
    return FALSE;
  }
  if(mod1->max_distance != mod2->max_distance ){
    return FALSE;
  }
  if(mod1->symbol != mod2->symbol ){
    return FALSE;
  }
  if(mod1->identifier != mod2->identifier ){
    return FALSE;
  }
  // aa list
  int i = 0;
  for(i=0; i<AA_LIST_LENGTH; i++){
    if(mod1->aa_list[i] != mod2->aa_list[i]){
      return FALSE;
    }
  }
  // everything matched
  return TRUE;
}


// FIXME: implement this
BOOLEAN_T is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod){

  MODIFIED_AA_T id = mod->identifier;
  if( (aa & id) != 0 ){  // should == id
    return TRUE;
  } 

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
  if( mod->aa_list[ (int)modified_aa_to_char(aa) - 'A' ] == TRUE ){
  //if( mod->aa_list[ (int)modified_aa_to_char(aa) ] == TRUE 
    //      && ! is_aa_modified(aa, mod)){
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
    carp(CARP_ERROR, "Cannot modify aa.  Either aa or mod is NULL.");
    return;
  }
  *aa = *aa | mod->identifier;
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

/**
 * \brief Generates a string representation of an aa_mod and returns a
 * pointer to that newly allocated string.
 */
char* aa_mod_to_string(AA_MOD_T* mod){
  const char* format_str = 
    "mass change=%.2f, symbol=%c, max=%d, position=%s, apply to ";
  char* return_str = (char*)mycalloc( strlen(format_str) + 50, sizeof(char));
  // add 26 for letters and some more for good measure

  // get position info
  const char* pos_format = "%c-term most %d from end";
  char* pos_buffer = (char*)mycalloc(strlen(pos_format), sizeof(char));
  switch(mod->position){
  case ANY_POSITION:
    strcpy(pos_buffer, "any");
    break;
  case C_TERM:
    sprintf(pos_buffer, pos_format, 'C', mod->max_distance);
    break;
  case N_TERM:
    sprintf(pos_buffer, pos_format, 'N', mod->max_distance);
    break;
  }
  int length = sprintf(return_str, format_str, mod->mass_change, 
                       mod->symbol, mod->max_per_peptide, pos_buffer);

  char* str_ptr = return_str + length;
  int i = 0;
  for(i=0; i<AA_LIST_LENGTH; i++){
    if( mod->aa_list[i] == TRUE ){
      sprintf(str_ptr, "%c", (i + 'A'));
      str_ptr++;
    }
  }
  free(pos_buffer);
  return return_str;
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

/**
 * \brief Create a string containing all of the amino acids that can
 * be modified by this aa_mod.  E.g. if S, T, and Y can be modified,
 * returns "STY".
 * \returns A newly allocated string.
 */
char* aa_mod_get_aa_list_string(AA_MOD_T* mod){

  char* list_string = (char*)mycalloc(AA_LIST_LENGTH+1, sizeof(char));
  int mod_idx, string_idx=0;
  for(mod_idx = 0; mod_idx < AA_LIST_LENGTH; mod_idx++){
    if(mod->aa_list[mod_idx] == TRUE){
      list_string[string_idx] = (char)('A' + mod_idx);
      string_idx++;
    }
    
  }// next aa
  list_string[string_idx] = '\0';
  return list_string;
}

/**
 * Count the number of modified aas in the string.
 */
int count_modified_aas(MODIFIED_AA_T* seq){
  if( seq == NULL ){
    return 0;
  }

  int count = 0;
  int aa_idx = 0;
  while(seq[aa_idx] != MOD_SEQ_NULL){
    if( GET_MOD_MASK & seq[aa_idx] ){
      count++;
    }
    aa_idx++;
  }

  return count;

}
