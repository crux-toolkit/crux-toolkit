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
#include "mass.h"
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
   MODIFIED_AA_T* and the index of an AA_MOD_T and return true if the
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
  bool aa_list[AA_LIST_LENGTH];
                       ///< an array indexed by AA, true if can be modified
  int max_per_peptide; ///< the maximum number of mods per peptide
  MOD_POSITION_T position; ///< where the mod can occur in the pep/prot
  int max_distance;        ///< the max distance from the protein terminus
  char symbol;         ///< the character to represent the mod in sqt files
  bool prevents_cleavage; ///< can modification prevent cleavage?
  bool prevents_xlink; ///< can modification prevent xlink?
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
  mod->prevents_cleavage = false;
  mod->prevents_xlink = false;
  // allocate the aa lists for mods 
  /*
  mod->aa_list =                        // all 0's?
    (bool*)mycalloc(AA_LIST_LENGTH, sizeof(bool)); 
  */
  // initialize to false
  int aa_idx = 0;
  for(aa_idx=0; aa_idx < AA_LIST_LENGTH; aa_idx++){
    mod->aa_list[aa_idx] = false;
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
  MODIFIED_AA_T mod_aa = (MODIFIED_AA_T)aa - (MODIFIED_AA_T)'A';
  return mod_aa;
}

/**
 * \brief Converts a MODIFIED_AA_T to it's textual representation,
 * i.e. a letter followed by between 0 and 11 symbols for the
 * modifications made to the amino acid.
 * \returns A newly allocated char* with amino acid and modifciation
 * symbols. 
 */
char* modified_aa_to_string_with_symbols(MODIFIED_AA_T aa){

  int modified_by = 0;
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for(int mod_idx = 0; mod_idx< total_mods; mod_idx++){
    if( is_aa_modified(aa, mod_list[mod_idx])){
      modified_by++;
    }
  }
  char* return_string = (char*)mymalloc((modified_by+2)*sizeof(char));
                                                  // 2 = aa + null term
  int return_idx = 0;
  return_string[return_idx] = modified_aa_to_char(aa);
  return_idx++;
  for(int mod_idx = 0; mod_idx < total_mods; mod_idx++){
    if( is_aa_modified(aa, mod_list[mod_idx])){
      return_string[return_idx] = aa_mod_get_symbol(mod_list[mod_idx]);
      return_idx++;
    }
  }
  return_string[return_idx] = '\0';
  return return_string;
}

/**
 * \brief Converts a MODIFIED_AA_T to it's textual representation,
 * i.e. a letter either alone or followed by square brackets containing
 * the mass(es) of any modifications. If mass_format is
 * MOD_MASSES_SEPARATE, all masses are listed in a comma-separated
 * list.  If MOD_MASS_ONLY, they are summed and returned in one
 * number.  If AA_PLUS_MOD, the mass of the residue plus the mass of
 * the modifciation(s) is printed. 
 * 
 * \returns A newly allocated char* with amino acid and modifciation
 * masses in square brackets.
 */
char* modified_aa_to_string_with_masses(MODIFIED_AA_T aa, 
                                        MASS_FORMAT_T mass_format,
                                        int precision){
  int modified_by = 0;
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for(int mod_idx = 0; mod_idx< total_mods; mod_idx++){
    if( is_aa_modified(aa, mod_list[mod_idx])){
      modified_by++;
    }
  }

  // return the char if no mods
  if( modified_by == 0 ){
    char* aa_str = (char*)mymalloc(2 * sizeof(char));
    sprintf(aa_str, "%c", modified_aa_to_char(aa));
    return aa_str;
  }

  // go through all mods and sum mass or append to string

  // at most 11 mods, each at most 7 char "000.00,", total 77
  char* mass_string = (char*)mymalloc(128 * sizeof(char));
  char* mass_string_ptr = mass_string;
  FLOAT_T summed_masses = 0;
  if( mass_format == AA_PLUS_MOD ){
    summed_masses = 
      get_mass_mod_amino_acid(aa, get_mass_type_parameter("isotopic-mass"));
  } else {

    for(int mod_idx = 0; mod_idx < total_mods; mod_idx++){
      if( is_aa_modified(aa, mod_list[mod_idx])){
        // either add to the sum or to the string
        if( mass_format == MOD_MASS_ONLY ){
          summed_masses += aa_mod_get_mass_change(mod_list[mod_idx]);
        } else {
          sprintf(mass_string_ptr, "%.*f,", precision,
                  aa_mod_get_mass_change(mod_list[mod_idx]));
          mass_string_ptr += strlen(mass_string_ptr);
        }
      }
    }
  }

  // combine aa and masses
  char* return_string = NULL;
  if( mass_format != MOD_MASSES_SEPARATE ){ // X[000.00]'/0'
    return_string = (char*)mymalloc(10 * sizeof(char));
    sprintf(return_string, "%c[%1.*f]", modified_aa_to_char(aa), 
            precision, summed_masses );
  } else { // X[000.00,etc]'/0'

    return_string = (char*)mymalloc((4 + strlen(mass_string)) * sizeof(char));
    // take off the last comma
    mass_string[(strlen(mass_string) - 1)] = '\0';
    sprintf(return_string, "%c[%s]", modified_aa_to_char(aa), mass_string);

  }
  free(mass_string);
  return return_string;
}

/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the mass change of
 * any modifications in brackets following the modified residue.  
 * If mass_format is MOD_MASSES_SEPARATE, all masses are listed in a
 * comma-separated list.  If MOD_MASS_ONLY, they are summed and
 * returned in one number.  If AA_PLUS_MOD, the mass of the residue
 * plus the mass of the modifciation(s) is printed.  
 *
 * \returns A newly allocated array of characters, a text
 * representation of the modified sequence.
 */
char* modified_aa_string_to_string_with_masses(
 MODIFIED_AA_T* aa_string, // the modified aa's to translate
 int length, // length of aa_string
 MASS_FORMAT_T mass_format) // which mass value to print
{
  if( aa_string == NULL || length == 0 ){
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }
  // get access to the mods
  AA_MOD_T** global_mod_list = NULL;
  int mod_list_length = get_all_aa_mod_list(&global_mod_list);

  // count up  modifications
  int count = 0;
  for(int mod_str_idx=0; mod_str_idx< length; mod_str_idx++){
    // check each mod in list
    for(int global_idx = 0; global_idx < mod_list_length; global_idx++){
      if( is_aa_modified(aa_string[mod_str_idx], 
                         global_mod_list[global_idx]) ){
        count++;  // count all mods
      }
    }
  }

  int precision = get_int_parameter("mod-precision");
  // max total length = #aas + ( #mods * (strlen("[000.,]")+precision) ) + '/0'
  int buffer_size = length + (count * (9 + precision)) + 1;
  char* return_string = (char*)mymalloc(buffer_size * sizeof(char));
  char* return_str_ptr = return_string;
  for(int mod_str_idx = 0; mod_str_idx<length; mod_str_idx++){

    char* cur_mod = 
      modified_aa_to_string_with_masses(aa_string[mod_str_idx], 
                                        mass_format, precision );
    strcpy( return_str_ptr, cur_mod );
    return_str_ptr += strlen(cur_mod);
    free(cur_mod);
  }

  return return_string;
}

/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the symbol for all
 * applied modifications.
 *
 * \returns A newly allocated array of characters, a text
 * representation of the modified sequence.
 */
char* modified_aa_string_to_string_with_symbols(
  MODIFIED_AA_T* aa_string, 
  int length)
{
  if( aa_string == NULL || length == 0 ){
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

    char* cur_mod = modified_aa_to_string_with_symbols(aa_string[mod_str_idx]);
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
 * \brief Allocates an array of MODIFIED_AA_T's and populates it with
 * the MODIFIED_AA_T value that corresponds to each sequence char
 * value and trailing modification symbols or masses.  Returns the new
 * sequence via the mod_sequence argument.
 *
 * \returns The length of the mod_sequence array.
 */
int convert_to_mod_aa_seq(const char* sequence, 
                          MODIFIED_AA_T** mod_sequence,
                          MASS_FORMAT_T mass_format){

  if( sequence == NULL ){
    carp(CARP_ERROR, "Cannot convert NULL sequence to modifiable characters"); 
    return 0;
  }

  MASS_TYPE_T mass_type = get_mass_type_parameter("isotopic-mass");

  int seq_len = strlen(sequence);
  MODIFIED_AA_T* new_sequence = 
    (MODIFIED_AA_T*)mycalloc( seq_len + 1, sizeof(MODIFIED_AA_T) );

  unsigned int seq_idx = 0;  // current position in given sequence
  unsigned int mod_idx = 0;  // current position in the new_sequence
  for(seq_idx = 0; seq_idx < strlen(sequence); seq_idx++){
    // is the character a residue?
    if( sequence[seq_idx] >= 'A' && sequence[seq_idx] <= 'Z' ){ 
      // add to the new sequence
      new_sequence[mod_idx] = char_aa_to_modified( sequence[seq_idx] );
      mod_idx++;
      continue;
    } 

    // else it's a modification as symbol or mass
    const AA_MOD_T* aa_mod = NULL;

    if( sequence[seq_idx] == '[' || sequence[seq_idx] == ','){//mod mass
      seq_idx++;
      FLOAT_T delta_mass = atof(sequence + seq_idx);
      if( mass_format == AA_PLUS_MOD ){
        assert(mod_idx > 0);
        delta_mass -= get_mass_mod_amino_acid(new_sequence[mod_idx - 1], 
                                              mass_type);
        }
      // translate mass into aa_mod
      aa_mod = get_aa_mod_from_mass(delta_mass);

      // move to next mod or residue
      while( sequence[seq_idx] != ']' && sequence[seq_idx] != ','){
        seq_idx++;
      }
      // back up one for the comma
      if( sequence[seq_idx] == ',' ){
        seq_idx--;
      }
    } else { // mod symbol
      // translate character into aa_mod
      aa_mod = get_aa_mod_from_symbol(sequence[seq_idx]);
    }
    if (mod_idx == 0) {
      //This can happen with nterminal modifications from comet.
      seq_idx++;
      if (seq_idx < strlen(sequence) && sequence[seq_idx] >= 'A'  && sequence[seq_idx] <= 'Z') {
        new_sequence[mod_idx] = char_aa_to_modified( sequence[seq_idx] );
        mod_idx++;
      } else {
        carp(CARP_FATAL, "Cannot parse sequence %s", sequence);
      }
      
      
    }
    // apply the modification
    if( aa_mod == NULL ){
      carp(CARP_WARNING, "There is an unidentifiable modification in sequence "
           "<%s> at position %d.", sequence, seq_idx - 1);
    } else {
      // apply modification
      modify_aa(&new_sequence[mod_idx-1], aa_mod);
    }
  } // next character in given sequence

  // null terminate
  new_sequence[mod_idx] = MOD_SEQ_NULL;

  //  return new_sequence;
  *mod_sequence = new_sequence;
  return mod_idx;
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
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(
  MODIFIED_AA_T* source ///< Sequence to copy
  ) {
  if (source == NULL) {
    carp(CARP_ERROR, "Cannot copy NULL sequence of modified_aa's.");
    return NULL;
  }

  size_t length = 0;
  while (source[length] != MOD_SEQ_NULL) {
    length++;
  }
  return copy_mod_aa_seq(source, length);

}


/**
 * \returns whether the two modified sequences are equal or not
 */
bool equal_seq(
  const MODIFIED_AA_T* seq1, ///< Sequence 1
  const MODIFIED_AA_T* seq2  ///< Sequence 2
  ) {

  size_t idx = 0;

  while (seq1[idx] != MOD_SEQ_NULL && seq2[idx] != MOD_SEQ_NULL) {
    if (seq1[idx] != seq2[idx]) {
      return false;
    }
    idx++;
  }

  return (seq1[idx] == MOD_SEQ_NULL && seq2[idx] == MOD_SEQ_NULL);

}

/**
 * \brief Remove any characters not A-Z from a peptide sequence.
 * \returns A newly allocated string with the given sequence less any
 * modififcation symbols or masses.
 */
char* unmodify_sequence(const char* modified_sequence){
  if( modified_sequence == NULL ){
    return NULL;
  }

  char* new_seq = my_copy_string(modified_sequence);
  unmodify_sequence_in_place(new_seq);
  return new_seq;

}

/**
 * \brief Remove any characters that are not A-Z from a given peptide
 * seqeunce.
 */
void unmodify_sequence_in_place(char* sequence){

  int seq_len = strlen(sequence);
  int next_aa_idx = 0;
  int move_to_here = 0;
  while( next_aa_idx < seq_len ){
    if( sequence[next_aa_idx] >= 'A' && sequence[next_aa_idx] <= 'Z' ){
      sequence[move_to_here] = sequence[next_aa_idx];
      move_to_here++;
    }
    next_aa_idx++;
  }
  sequence[move_to_here] = '\0';
}


/**
 * \brief Determine if an array of MODIFIED_AA_T is a palindrome.  
 * Used by reverse_sequence to avoid returning a reversed sequence
 * that is the same as the target.  Ignores the first and last
 * residues. 
 * \returns true if the reversed sequence would be the same as the
 * forward, otherwise false.
 */
bool modified_aa_seq_is_palindrome(MODIFIED_AA_T* seq, int length){
  if( seq == NULL ){
    return false;
  }

  int left_idx = 1;    // skip first and last residues
  int right_idx = length - 2;

  while( left_idx < right_idx ){
    if( seq[left_idx] != seq[right_idx] ){
      return false;
    }// else, keep checking
    left_idx++;
    right_idx--;
  }

  // if we got to here, they all matched
  return true;
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
 * return false.  If the given mod list does not exactly match the
 * mods read from the parameter file (including the order in which
 * they are listed) return false.  If returning false, print a warning
 * with the lines that should be included in the parameter file.
 *
 * \returns true if the given mods are the same as those from the
 * parameter file.
 */
bool compare_mods(AA_MOD_T** psm_file_mod_list, int file_num_mods){
  AA_MOD_T** mod_list = NULL;
  int num_mods = get_all_aa_mod_list(&mod_list);

  if( num_mods != file_num_mods ){
    return false;
  }
  int mod_idx = 0;
  for(mod_idx=0; mod_idx<num_mods; mod_idx++){
    if( ! compare_two_mods(mod_list[mod_idx], psm_file_mod_list[mod_idx]) ){
      print_a_mod(mod_list[mod_idx]);
      print_a_mod(psm_file_mod_list[mod_idx]);
      return false;
    }
  }
  return true;
}

/**
 * \brief Compare two mods to see if they are the same, i.e. same mass
 * change, unique identifier, position
 */
bool compare_two_mods(AA_MOD_T* mod1, AA_MOD_T* mod2){
  if(mod1->mass_change != mod2->mass_change ){
    return false;
  }
  if(mod1->max_per_peptide != mod2->max_per_peptide ){
    return false;
  }
  if(mod1->position != mod2->position ){
    return false;
  }
  if(mod1->max_distance != mod2->max_distance ){
    return false;
  }
  if(mod1->symbol != mod2->symbol ){
    return false;
  }
  if(mod1->identifier != mod2->identifier ){
    return false;
  }
  // aa list
  int i = 0;
  for(i=0; i<AA_LIST_LENGTH; i++){
    if(mod1->aa_list[i] != mod2->aa_list[i]){
      return false;
    }
  }
  // everything matched
  return true;
}


// FIXME: implement this
bool is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod){

  MODIFIED_AA_T id = mod->identifier;
  if( (aa & id) != 0 ){  // should == id
    return true;
  } 

  return false;
}

/**
 * \brief Determine if this modified amino acid can be modified by
 * this modification.
 *
 * Checks the mod list of modifiable residues to see if aa is in the
 * list.  Also checks to see if aa has already been modified by this
 * mod.  
 * \returns true if it can be modified, else false
 */
bool is_aa_modifiable
 (MODIFIED_AA_T aa, ///< The sequence amino acid to be modified
  AA_MOD_T* mod){   ///< The aa_mod to modify it with

  if( is_aa_modified(aa, mod) == true ){
    return false;
  }
  if( mod->aa_list[ (int)modified_aa_to_char(aa) - 'A' ] == true ){
  //if( mod->aa_list[ (int)modified_aa_to_char(aa) ] == true 
    //      && ! is_aa_modified(aa, mod)){
    return true;
  }
  // else not in list or already modified by this mod
  return false;
}

/**
 * \brief Adds a modification to a MODIFIED_AA_T.
 *
 * Assumes that the aa is modifiable, no explicit check.  If the aa is
 * already modified for the mod, no change to aa.
 */
void modify_aa(MODIFIED_AA_T* aa, const AA_MOD_T* mod){
  if( aa == NULL || mod == NULL ){
    carp(CARP_ERROR, "Cannot modify aa.  Either aa or mod is NULL.");
    return;
  }
  *aa = *aa | mod->identifier;
}

/**
 * \brief Return the AA_MOD_T associated with the given symbol.  If
 * the symbol does not represent a modification, returns null.
 * Requires that parameters have been initialized.
 */
const AA_MOD_T* get_aa_mod_from_symbol(const char symbol){

  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for(int mod_idx = 0; mod_idx < total_mods; mod_idx++){
    AA_MOD_T* cur_mod = mod_list[mod_idx];
    if( cur_mod->symbol == symbol ){
      return cur_mod;
    }
  }

  // none of the mods matched the symbol
  carp(CARP_ERROR, "The modification symbol '%c' is not valid.", symbol);
  return NULL;
}

/**
 * \brief Return the delta mass associated with the given modification
 * symbol.  If the symbol does not represent a modification, returns
 * 0. Requires that parameters have been initialized.
 */
FLOAT_T get_mod_mass_from_symbol(const char symbol){

  const AA_MOD_T* mod = get_aa_mod_from_symbol(symbol);
  if( mod == NULL ){
    return 0;
  }

  return aa_mod_get_mass_change(mod);
}

// a temporary AA_MOD_T with identifier and mass for (potentially)
// more than one modification
AA_MOD_T multi_mod; 

/**
 * \brief Return the modification identifier that can be used to modify a
 * MODIFIED_AA_T so that it has the given mass shift.
 * The mass may either be from a single AA_MOD_T as given by the user
 * or from any combination of AA_MOD_T's.  If no combinations of
 * AA_MOD_T's can be found for the mass, returns 0.
 * Requires that parameters have been initialized.
 */
const AA_MOD_T* get_aa_mod_from_mass(FLOAT_T mass){

  // find the identifier for this mass shift
  MODIFIED_AA_T id = get_mod_identifier(mass);

  if( id == 0 ){ // get_mod_id already warned
    //here we are creating a modification.  Hopefully
    //this will only get used when post-processing search
    //results, maybe we will need to check that at some point
    //SJM 2013_10_10.
    carp(CARP_INFO, "Creating modification for %f", mass);
    AA_MOD_T** mods = NULL;
    int num_mods = get_all_aa_mod_list(&mods);
    mods[num_mods]->mass_change = mass;
    incrementNumMods(); 
    initialize_aa_mod_combinations_array();
    return get_aa_mod_from_mass(mass);
  } else {
    // set multi_mod_identifier to that id and mass
    multi_mod.identifier = id;
    
    
    
  }
  multi_mod.mass_change = mass;
  carp(CARP_DETAILED_DEBUG, "Returning %i %f", multi_mod.identifier, multi_mod.mass_change);
  return &multi_mod;
}

/**
 * print all fields in aa mod.  For debugging
 */
void print_a_mod(AA_MOD_T* mod){
  printf("AMOD: mass %.2f, max per %d, max dist %d, symb %c, aa list ",
         mod->mass_change, mod->max_per_peptide, mod->max_distance,
         mod->symbol);

  int i = 0;
  for(i=0; i<AA_LIST_LENGTH; i++){
    if( mod->aa_list[i] == true ){
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
  char* pos_buffer = (char*)mycalloc(strlen(pos_format) + 9, sizeof(char));
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
    if( mod->aa_list[i] == true ){
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
double aa_mod_get_mass_change(const AA_MOD_T* mod){
  return mod->mass_change;
}

/**
 * \brief Access to the aa_list of the AA_MOD_T struct.  This pointer
 * can be used to get or set the list of residues on which this mod
 * can be placed.
 * \returns A pointer to the list of amino acids on which this mod can
 * be placed.
 */
bool* aa_mod_get_aa_list(AA_MOD_T* mod){
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
 * \brief Sets whether the modification can prevent cleavage.
 * \returns void
 */
void aa_mod_set_prevents_cleavage(AA_MOD_T* mod, bool prevents_cleavage) {
  mod->prevents_cleavage=prevents_cleavage;
}

/**
 * \brief gets whether the modification can prevent cleavage
 * \returns true or false
 */
bool aa_mod_get_prevents_cleavage(AA_MOD_T* mod) {
  return mod->prevents_cleavage;
}

/**
 * \brief Sets whether the modifications can prevent cross-linking.
 * \returns void
 */
void aa_mod_set_prevents_xlink(AA_MOD_T* mod, bool prevents_xlink) {
  mod->prevents_xlink = prevents_xlink;
}

/**
 * \brief gets whether the modification can prevent cross-linking.
 * \returns true or false
 */
bool aa_mod_get_prevents_xlink(AA_MOD_T* mod) {
  return mod->prevents_xlink;
}

/**
 * \brief The character used to uniquely identify the mod in the sqt file.
 * \returns The character identifier.
 */
char aa_mod_get_symbol(const AA_MOD_T* mod){
  return mod->symbol;
}

/**
 * \brief The bitmask used to uniquely identify the mod.
 * \returns The short int bitmask used to identify the mod.
 */
int aa_mod_get_identifier(const AA_MOD_T* mod){
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
    if(mod->aa_list[mod_idx] == true){
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

/**
 * /returns the mass of the modified sequence
 */
FLOAT_T get_mod_aa_seq_mass(
  MODIFIED_AA_T* seq, ///< The modified sequence
  MASS_TYPE_T mass_type ///<  mono or average?
  ) {
  
  // get access to the mods
  AA_MOD_T** global_mod_list = NULL;
  int mod_list_length = get_all_aa_mod_list(&global_mod_list);
  
  if (seq == NULL) {
    return 0.0;
  }
  
  FLOAT_T ans = 0;
  int aa_idx = 0;
  while(seq[aa_idx] != MOD_SEQ_NULL) {
    ans += get_mass_amino_acid(modified_aa_to_char(seq[aa_idx]), mass_type);
    if ( GET_MOD_MASK & seq[aa_idx]) {
      for(int global_idx = 0; global_idx < mod_list_length; global_idx++){
        if( is_aa_modified(seq[aa_idx], 
                     global_mod_list[global_idx]) ){
          ans += aa_mod_get_mass_change(global_mod_list[global_idx]);
        }
      }
    }
    aa_idx++;
  }
  if (mass_type == AVERAGE) {
    return ans + MASS_H2O_AVERAGE;
  } else {
    return ans + MASS_H2O_MONO;
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
