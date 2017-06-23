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

#include "model/Peptide.h"
#include "modifications.h"
#include "mass.h"
#include "GlobalParams.h"
#include "Params.h"
#include <stack>
#include <iostream>

using namespace std;

/**
 * \brief Keeps a stack of allocated modification sequences
 */
class MODIFIED_AA_T_Cache {
 protected:
  stack<MODIFIED_AA_T*> cache_; ///< Cache itself
 public:
  
  /**
   * Constructor
   */
  MODIFIED_AA_T_Cache() {;}

  /**
   * Destructor, deletes all objects in cache
   */
  ~MODIFIED_AA_T_Cache() {
    size_t objects_in_cache = cache_.size();
    while(!cache_.empty()) {
      MODIFIED_AA_T *element = cache_.top();
      cache_.pop();
      delete []element;
    }
  }

  /**
   * checks out a object out of the cache
   */
  MODIFIED_AA_T* checkout(bool clear=false) {
    MODIFIED_AA_T* new_element;
    if (cache_.empty()) {
      new_element = new MODIFIED_AA_T[MAX_PEPTIDE_LENGTH + 1];
    } else {
      new_element = cache_.top();
      cache_.pop();
    }
    if (clear) {
      memset(new_element, 0, sizeof(MODIFIED_AA_T)*(MAX_PEPTIDE_LENGTH+1));
    }
    return (new_element);
  }
  
  
  /**
   * checks in a object into the cache
   */
  void checkin(
    MODIFIED_AA_T* array
  ) {
    cache_.push(array);
  }
  
};

MODIFIED_AA_T_Cache modified_aa_cache; ///< The cache of modified sequences

/**
 * \return a new modification sequence array. Can be off of already allocated
 * arrays
 */
MODIFIED_AA_T* newModSeq() {
  MODIFIED_AA_T* ans = modified_aa_cache.checkout();
  return(ans);
}

/**
 * Releases the modification sequence back into the cache
 */
void freeModSeq(
  MODIFIED_AA_T* &seq  ///< sequence to release
  ) {
  modified_aa_cache.checkin(seq);
  seq=NULL;
}

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

uint16_t mod_id_masks[MAX_AA_MODS] = {
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

/* Definitions of public methods */

/**
 * \brief Allocate an AA_MOD, including the aa_list and initialize all
 * fields to default values.  Symbol and unique identifier are set
 * according to index.  
 * \returns A heap allocated AA_MOD with default values.
 */
AA_MOD_T::AA_MOD_T(int mod_idx) {
  assert(mod_idx > -1 && mod_idx < MAX_AA_MODS );

  mass_change_ = 0;
  max_per_peptide_ = 0;
  position_ = ANY_POSITION;
  max_distance_ = MAX_PROTEIN_SEQ_LENGTH;
  symbol_ = mod_sqt_symbols[mod_idx];
  identifier_ = mod_id_masks[mod_idx];
  prevents_cleavage_ = false;
  prevents_xlink_ = false;
  mono_link_ = false;
  // initialize to false
  for (int aa_idx = 0; aa_idx < AA_LIST_LENGTH; aa_idx++) {
    aa_list_[aa_idx] = false;
  }
}

/**
 * \brief Free the memory for an AA_MOD including the aa_list.
 */
AA_MOD_T::~AA_MOD_T() {
  ;
}

/**
 * Set bitmask id
 */
void AA_MOD_T::setModIdx(int mod_idx) {
  assert(mod_idx > -1 && mod_idx < MAX_AA_MODS );
  identifier_ = mod_id_masks[mod_idx];
  symbol_ = mod_sqt_symbols[mod_idx];
}
  
/**
 * \returns bitmask id
 */
int AA_MOD_T::getIdentifier() {
  return(identifier_);
}

void AA_MOD_T::setIdentifier(int identifier) {
  identifier_ = identifier;
}

/**
 * \brief Converts a MODIFIED_AA into a char, effectively unmodifying it.
 * \returns The unmodified char representation of an aa.
 */
char modified_aa_to_char(MODIFIED_AA_T aa) {
  aa = aa & 0x001F;         // 0000 0000 0001 1111
  return (char)aa + 'A';
}

/**
 * \brief Converts a char representation of an aa to an unmodified aa
 * of type MODIFIED_AA_T.  Requires 'A' <= char <= 'Z'. 
 * \returns The MODIFIED_AA_T represnetation of an aa.
 */
MODIFIED_AA_T char_aa_to_modified(char aa) {
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
char* modified_aa_to_string_with_symbols(MODIFIED_AA_T aa) {

  int modified_by = 0;
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < total_mods; mod_idx++) {
    if (mod_list[mod_idx]->isModified(aa)) {
      modified_by++;
    }
  }
  char* return_string = (char*)mymalloc((modified_by+2)*sizeof(char));
                                                  // 2 = aa + null term
  int return_idx = 0;
  return_string[return_idx] = modified_aa_to_char(aa);
  return_idx++;
  for (int mod_idx = 0; mod_idx < total_mods; mod_idx++) {
    if (mod_list[mod_idx]->isModified(aa)) {
      return_string[return_idx] = mod_list[mod_idx]->getSymbol();
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
                                        int precision) {
  int modified_by = 0;
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < total_mods; mod_idx++) {
    if (mod_list[mod_idx]->isModified(aa)) {
      modified_by++;
    }
  }

  // return the char if no mods
  if (modified_by == 0) {
    char* aa_str = (char*)mymalloc(2 * sizeof(char));
    sprintf(aa_str, "%c", modified_aa_to_char(aa));
    return aa_str;
  }

  // go through all mods and sum mass or append to string

  // at most 11 mods, each at most 7 char "000.00,", total 77
  char* mass_string = (char*)mymalloc(128 * sizeof(char));
  char* mass_string_ptr = mass_string;
  FLOAT_T summed_masses = 0;
  if (mass_format == AA_PLUS_MOD) {
    summed_masses = 
      get_mass_mod_amino_acid(aa, GlobalParams::getIsotopicMass());
  } else {

    for (int mod_idx = 0; mod_idx < total_mods; mod_idx++) {
      if (mod_list[mod_idx]->isModified(aa)) {
        // either add to the sum or to the string
        if (mass_format == MOD_MASS_ONLY) {
          summed_masses += mod_list[mod_idx]->getMassChange();
        } else {
          sprintf(mass_string_ptr, "%.*f,", precision,
                  mod_list[mod_idx]->getMassChange());
          mass_string_ptr += strlen(mass_string_ptr);
        }
      }
    }
  }

  // combine aa and masses
  char* return_string = NULL;
  if (mass_format != MOD_MASSES_SEPARATE) { // X[000.00]'/0'
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
  MASS_FORMAT_T mass_format // which mass value to print
) {
  if (aa_string == NULL || length == 0) {
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }
  // get access to the mods
  AA_MOD_T** global_mod_list = NULL;
  int mod_list_length = get_all_aa_mod_list(&global_mod_list);

  // count up  modifications
  int count = 0;
  for (int mod_str_idx = 0; mod_str_idx < length; mod_str_idx++) {
    // check each mod in list
    for(int global_idx = 0; global_idx < mod_list_length; global_idx++) {
      if (global_mod_list[global_idx]->isModified(aa_string[mod_str_idx])) {
        count++;  // count all mods
      }
    }
  }

  int precision = GlobalParams::getModPrecision();
  // max total length = #aas + ( #mods * (strlen("[000.,]")+precision) ) + '/0'
  int buffer_size = length + (count * (9 + precision)) + 1;
  char* return_string = (char*)mymalloc(buffer_size * sizeof(char));
  char* return_str_ptr = return_string;
  for (int mod_str_idx = 0; mod_str_idx < length; mod_str_idx++) {
    char* cur_mod = 
      modified_aa_to_string_with_masses(aa_string[mod_str_idx], mass_format, precision);
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
  int length
) {
  if (aa_string == NULL || length == 0) {
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }

  AA_MOD_T** global_mod_list = NULL;
  int mod_list_length = get_all_aa_mod_list(&global_mod_list);
  int global_idx = 0;

  // count up amino acids and modifications
  int count = 0;
  int mod_str_idx = 0;
  //  while (aa_string[mod_str_idx] != MOD_SEQ_NULL) {
  for (mod_str_idx = 0; mod_str_idx < length; mod_str_idx++) {
    count++;      // count the aas
    for (global_idx = 0; global_idx < mod_list_length; global_idx++) {
      if (global_mod_list[global_idx]->isModified(aa_string[mod_str_idx])) {
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
  //  for (mod_str_idx = 0; mod_str_idx < mod_str_len; mod_str_idx++) {
  for (mod_str_idx = 0; mod_str_idx < length; mod_str_idx++) {

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
char* modified_aa_to_unmodified_string(MODIFIED_AA_T* aa_string, int length) {
  if (aa_string == NULL) {
    carp(CARP_ERROR, "Cannot print a NULL modified sequence");
    return NULL;
  }

  char* new_string = (char*)mycalloc(length+1, sizeof(char));
  for (int aa_idx = 0; aa_idx < length; aa_idx++) {
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
int convert_to_mod_aa_seq(const string& sequence, 
                          MODIFIED_AA_T** mod_sequence,
                          MASS_FORMAT_T mass_format) {

  if( sequence.empty() ){
    carp(CARP_ERROR, "Cannot convert NULL sequence to modifiable characters"); 
    return 0;
  }

  MASS_TYPE_T mass_type = GlobalParams::getIsotopicMass();
  const char* csequence = sequence.c_str();

  int seq_len = sequence.length();
  MODIFIED_AA_T* new_sequence = newModSeq(); 
  unsigned int seq_idx = 0;  // current position in given sequence
  unsigned int mod_idx = 0;  // current position in the new_sequence
  for(seq_idx = 0; seq_idx < seq_len; seq_idx++){
    // is the character a residue?
    if( csequence[seq_idx] >= 'A' && csequence[seq_idx] <= 'Z' ){ 
      // add to the new sequence
      new_sequence[mod_idx] = char_aa_to_modified( csequence[seq_idx] );
      mod_idx++;
      continue;
    } 

    // else it's a modification as symbol or mass
    const AA_MOD_T* aa_mod = NULL;

    if (sequence[seq_idx] == '[' || sequence[seq_idx] == ',') {//mod mass
      seq_idx++;
      FLOAT_T delta_mass = atof(csequence + seq_idx);
      if (mass_format == AA_PLUS_MOD) {
        assert(mod_idx > 0);
        delta_mass -= get_mass_mod_amino_acid(new_sequence[mod_idx - 1], 
                                              mass_type);
      }
      // translate mass into aa_mod
      aa_mod = get_aa_mod_from_mass(delta_mass);

      // move to next mod or residue
      while (sequence[seq_idx] != ']' && sequence[seq_idx] != ',') {
        seq_idx++;
      }
      // back up one for the comma
      if (sequence[seq_idx] == ',') {
        seq_idx--;
      }
    } else { // mod symbol
      // translate character into aa_mod
      aa_mod = get_aa_mod_from_symbol(sequence[seq_idx]);
    }
    if (mod_idx == 0) {
      //This can happen with nterminal modifications from comet.
      seq_idx++;
      if (seq_idx < seq_len && sequence[seq_idx] >= 'A'  && sequence[seq_idx] <= 'Z') {
        new_sequence[mod_idx] = char_aa_to_modified( sequence[seq_idx] );
        mod_idx++;
      } else {
        carp(CARP_FATAL, "Cannot parse sequence %s", csequence);
      }
      
      
    }
    // apply the modification
    if (aa_mod == NULL) {
      carp(CARP_WARNING, "There is an unidentifiable modification in sequence "
           "<%s> at position %d.", csequence, seq_idx - 1);
    } else {
      // apply modification
      aa_mod->modify(&new_sequence[mod_idx-1]);
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
MODIFIED_AA_T* copy_mod_aa_seq(const MODIFIED_AA_T* source, int length) {
  if (source == NULL) {
    carp(CARP_ERROR, "Cannot copy NULL sequence of modified_aa's.");
    return NULL;
  }

  MODIFIED_AA_T* new_seq = newModSeq();
  memcpy( new_seq, source, length * sizeof(MODIFIED_AA_T));
  new_seq[length] = MOD_SEQ_NULL;

  return new_seq;
}

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(
  const MODIFIED_AA_T* source ///< Sequence to copy
  ) {
  if (source == NULL) {
    carp(CARP_ERROR, "Cannot copy NULL sequence of modified_aa's.");
    return NULL;
  }

  MODIFIED_AA_T* new_seq = newModSeq();
  size_t length = 0;
  while (source[length] != MOD_SEQ_NULL) {
    new_seq[length] = source[length];
    length++;
  }
  new_seq[length] = MOD_SEQ_NULL;
  return(new_seq);

}

/**
 * \brief Determine if an array of MODIFIED_AA_T is a palindrome.  
 * Used by reverse_sequence to avoid returning a reversed sequence
 * that is the same as the target.  Ignores the first and last
 * residues. 
 * \returns true if the reversed sequence would be the same as the
 * forward, otherwise false.
 */
bool modified_aa_seq_is_palindrome(MODIFIED_AA_T* seq, int length) {
  if (seq == NULL) {
    return false;
  }

  int left_idx = 1;    // skip first and last residues
  int right_idx = length - 2;

  while (left_idx < right_idx) {
    if (seq[left_idx] != seq[right_idx]) {
      return false;
    }// else, keep checking
    left_idx++;
    right_idx--;
  }

  // if we got to here, they all matched
  return true;
}

// FIXME: implement this
bool AA_MOD_T::isModified(MODIFIED_AA_T aa) {
  return (aa & identifier_) != 0;
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
bool AA_MOD_T::isModifiable(
  MODIFIED_AA_T aa ///< The sequence amino acid to be modified
  ) { 

  if (isModified(aa)) {
    return false;
  }
  if (aa_list_[(int)modified_aa_to_char(aa) - 'A']) {
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
void AA_MOD_T::modify(MODIFIED_AA_T* aa) const {
  if (aa == NULL) {
    carp(CARP_ERROR, "Cannot modify aa.  aa is NULL.");
    return;
  }
  *aa = *aa | identifier_;
}

/**
 * \brief Return the AA_MOD_T associated with the given symbol.  If
 * the symbol does not represent a modification, returns null.
 * Requires that parameters have been initialized.
 */
const AA_MOD_T* get_aa_mod_from_symbol(const char symbol) {
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  for (int mod_idx = 0; mod_idx < total_mods; mod_idx++) {
    AA_MOD_T* cur_mod = mod_list[mod_idx];
    if (cur_mod->getSymbol() == symbol) {
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
FLOAT_T get_mod_mass_from_symbol(const char symbol) {

  const AA_MOD_T* mod = get_aa_mod_from_symbol(symbol);
  if (mod == NULL) {
    return 0;
  }

  return mod->getMassChange();
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
const AA_MOD_T* get_aa_mod_from_mass(FLOAT_T mass) {
  // find the identifier for this mass shift
  MODIFIED_AA_T id = get_mod_identifier(mass);

  if (id == 0) { // get_mod_id already warned
    //here we are creating a modification.  Hopefully
    //this will only get used when post-processing search
    //results, maybe we will need to check that at some point
    //SJM 2013_10_10.
    carp(CARP_INFO, "Creating modification for %f", mass);
    AA_MOD_T** mods = NULL;
    int num_mods = get_all_aa_mod_list(&mods);
    mods[num_mods]->setMassChange(mass);
    incrementNumMods(); 
    initialize_aa_mod_combinations_array();
    return get_aa_mod_from_mass(mass);
  } else {
    // set multi_mod_identifier to that id and mass
    multi_mod.setIdentifier(id);
  }
  multi_mod.setMassChange(mass);
  return &multi_mod;
}

/**
 * print all fields in aa mod.  For debugging
 */
void AA_MOD_T::print() {
  printf("AMOD: mass %.2f, max per %d, max dist %d, symb %c, aa list ",
         mass_change_, max_per_peptide_, max_distance_,
         symbol_);

  int i = 0;
  for (i = 0; i < AA_LIST_LENGTH; i++) {
    if (aa_list_[i]) {
      printf("%c", (i + 'A'));
    }
  }
  if (prevents_xlink_) {
    printf(" prevents-xlink");
  }
  if (prevents_cleavage_) {
    printf(" prevents-cleavage");
  }
  if (mono_link_) {
    printf(" mono");
  }
  printf("\n");
}

/**
 * \brief Generates a string representation of an aa_mod and returns a
 * pointer to that newly allocated string.
 */
char* AA_MOD_T::toCString() {
  const char* format_str = 
    "mass change=%.2f, symbol=%c, max=%d, position=%s, apply to ";
  char* return_str = (char*)mycalloc( strlen(format_str) + 50, sizeof(char));
  // add 26 for letters and some more for good measure

  // get position info
  const char* pos_format = "%c-term most %d from end";
  char* pos_buffer = (char*)mycalloc(strlen(pos_format) + 9, sizeof(char));
  switch (position_) {
  case ANY_POSITION:
    strcpy(pos_buffer, "any");
    break;
  case C_TERM:
    sprintf(pos_buffer, pos_format, 'C', max_distance_);
    break;
  case N_TERM:
    sprintf(pos_buffer, pos_format, 'N', max_distance_);
    break;
  }
  int length = sprintf(return_str, format_str, mass_change_, 
                       symbol_, max_per_peptide_, pos_buffer);

  char* str_ptr = return_str + length;
  int i = 0;
  for (i = 0; i < AA_LIST_LENGTH; i++) {
    if (aa_list_[i]) {
      sprintf(str_ptr, "%c", (i + 'A'));
      str_ptr++;
    }
  }
  free(pos_buffer);
  return return_str;
}

/* Setters and Getters */

/**
 * \brief Set the symbol used by this modification.
 * \returns void
 */
void AA_MOD_T::setSymbol(char symbol) {
  symbol_ = symbol;
}

/**
 * \brief Set the mass change caused by this modification.
 * \returns void
 */
void AA_MOD_T::setMassChange(double mass_change) {
  mass_change_ = mass_change;
}

/**
 * \brief Get the mass change caused by this modification.
 * \returns The mass change caused by this modification.
 */
double AA_MOD_T::getMassChange() const {
  return mass_change_;
}

/**
 * \brief Access to the aa_list of the AA_MOD_T struct.  This pointer
 * can be used to get or set the list of residues on which this mod
 * can be placed.
 * \returns A pointer to the list of amino acids on which this mod can
 * be placed.
 */
bool* AA_MOD_T::getAAList() {
  return aa_list_;
}

void AA_MOD_T::setAminoAcids(const char* aas, int num_aas) {
  for (int idx = 0;idx < num_aas;idx++) {
    setAminoAcid(aas[idx]);
  }
}

void AA_MOD_T::setAminoAcid(const char aa) {
  aa_list_[aa - 'A'] = true;
}


/**
 * \brief Set the maximum number of times this modification can be
 * placed on one peptide.
 * \returns void
 */
void AA_MOD_T::setMaxPerPeptide(int max) {
  max_per_peptide_ = max;
}

/**
 * \brief Get the maximum number of times this modification can be
 * placed on one peptide.  
 * \returns The max times per peptide this mod can be placed.
 */
int AA_MOD_T::getMaxPerPeptide() {
  return max_per_peptide_;
}

/**
 * \brief Set the maximum distance from the protein terminus that the
 * modification can be placed.  Which terminus (C or N) is determined
 * by the position type.  To indicate no position restriction, set to
 * MAX_PROTEIN_SEQ_LENGTH. 
 * \returns void
 */
void AA_MOD_T::setMaxDistance(int distance) {
  max_distance_ = (distance == -1) ? MAX_PROTEIN_SEQ_LENGTH : 
                                         distance;
}

/**
 * \brief Get the maximum distance from the protein end that the
 * modification can be placed.  Will be MAX_PROTEIN_SEQ_LENGTH if
 * position type is ANY_POSITION.
 * \returns Maximum distance from protein terminus at which mod can be
 * placed. 
 */
int AA_MOD_T::getMaxDistance() {
  return max_distance_;
}

/**
 * \brief Set the position type of an aa_mod.
 * \returns void
 */
void AA_MOD_T::setPosition(MOD_POSITION_T position) {
  position_ = position;
}

/**
 * \brief Where in the peptide can the modification be placed.
 * \returns ANY_POSITION for standard mods; C_TERM or N_TERM for those
 * that can only be placed at the ends of the peptide.
 */
MOD_POSITION_T AA_MOD_T::getPosition() {
  return position_;
}

/**
 * \brief Sets whether the modification can prevent cleavage.
 * \returns void
 */
void AA_MOD_T::setPreventsCleavage(bool prevents_cleavage) {
  prevents_cleavage_ = prevents_cleavage;
}

/**
 * \brief gets whether the modification can prevent cleavage
 * \returns true or false
 */
bool AA_MOD_T::getPreventsCleavage() {
  return prevents_cleavage_;
}

/**
 * \brief Sets whether the modifications can prevent cross-linking.
 * \returns void
 */
void AA_MOD_T::setPreventsXLink(bool prevents_xlink) {
  prevents_xlink_ = prevents_xlink;
}

/**
 * \brief gets whether the modification can prevent cross-linking.
 * \returns true or false
 */
bool AA_MOD_T::getPreventsXLink() {
  return prevents_xlink_;
}

/**
 * \brief tells whether the modification is a mono link.  Used by sfx
 */
void AA_MOD_T::setMonoLink(bool mono_link) {
  mono_link_ = mono_link;
}

bool AA_MOD_T::getMonoLink() {
  return(mono_link_);
}

/**
 * \brief The character used to uniquely identify the mod in the sqt file.
 * \returns The character identifier.
 */
char AA_MOD_T::getSymbol() {
  return symbol_;
}

/**
 * \brief Create a string containing all of the amino acids that can
 * be modified by this aa_mod.  E.g. if S, T, and Y can be modified,
 * returns "STY".
 * \returns A newly allocated string.
 */
string AA_MOD_T::getAAListString() {
  string s;
  for (int i = 0; i < AA_LIST_LENGTH; i++) {
    if (aa_list_[i]) {
      s += (char)('A' + i);
    }
  }
  return s;
}

/**
 * Count the number of modified aas in the string.
 */
int count_modified_aas(MODIFIED_AA_T* seq) {
  if (seq == NULL) {
    return 0;
  }

  int count = 0;
  int aa_idx = 0;
  while (seq[aa_idx] != MOD_SEQ_NULL) {
    if (GET_MOD_MASK & seq[aa_idx]) {
      count++;
    }
    aa_idx++;
  }

  return count;

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
