/**
 * \file modifications.h
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
 * $Revision: 1.7 $
 */
#ifndef MODIFICATION_FILE_H
#define MODIFICATION_FILE_H

#include <assert.h>
#include "utils.h"
#include "linked_list.h"
#include "objects.h"
#include "parameter.h"

/* Public constants */
static const int MAX_AA_MODS = 11;
static const int MAX_PROTEIN_SEQ_LENGTH = 40000;
static const int AA_LIST_LENGTH = 26; // A-Z
static const int MOD_MASS_PRECISION = 2;  // printed as X[1.23]X
#define GET_AA_MASK  0x001F   // 0000 0000 0001 1111
#define GET_MOD_MASK 0xFFE0   // 1111 1111 1110 0000


// this was moved to object.h b/c methods in peptide.h weren't compiling
//typedef unsigned short MODIFIED_AA_T; ///< letters in the expanded peptide
                                      ///alphabet, bits for mod1 mod2...aa
#define MOD_SEQ_NULL (MODIFIED_AA_T)('Z' - 'A' + 1) 
//enum {MOD_SEQ_NULL = 'Z' - 'A' + 1}; 
///< null terminating character of array

/*
   e.g. There are three aa_mods specified in this run and they are
   given the masks mod1  1000_0000_0000_0000
                   mod2  0100_0000_0000_0000
                   mod3  0010_0000_0000_0000
   The amino acid Y has the value
                      Y  0000_0000_0001_0011

   Suppose you set the variable aa to Y.  Do stuff.  To ask "is aa
  modified by mod2", do this. answer = mod2 AND aa; answer == 0
   To change aa to be modified by modcalle2, aa = mod2 NOR aa
 (i.e. !(mod2 || aa) )  now aa == 0100_0000_0001_0011
   If we ask again, answer == 0100_0000_0000_0000 
 */

/**
 * \brief Allocate an AA_MOD, including space for the aa_list and
 * initialize all fields to default values.  Symbol and unique
 * identifier are set according to index.  
 * \returns A heap allocated AA_MOD with default values.
 */
AA_MOD_T* new_aa_mod(int mod_idx);

/**
 * \brief Free the memory for an AA_MOD including the aa_list.
 */
void free_aa_mod(AA_MOD_T*);

/**
 * \brief Converts a MODIFIED_AA into a char, effectively unmodifying it.
 * \returns The unmodified char representation of an aa.
 */
char modified_aa_to_char(MODIFIED_AA_T aa);

/**
 * \brief Converts a char representation of an aa to an unmodified aa
 * of type MODIFIED_AA_T.  Requires 'A' <= char <= 'Z'. 
 * \returns The MODIFIED_AA_T represnetation of an aa.
 */
MODIFIED_AA_T char_aa_to_modified(char aa);

/**
 * \brief Converts a MODIFIED_AA_T* to it's textual representation,
 * i.e. a letter followed by between 0 and 11 symbols for the
 * modifications made to the amino acid.
 * \returns A newly allocated char* with amino acid and modifciation
 * symbols. 
 */
char* modified_aa_to_string_with_symbols(MODIFIED_AA_T aa);

/**
 * \brief Converts a MODIFIED_AA_T to it's textual representation,
 * i.e. a letter either alone or followed by square braces containing
 * the mass(es) of any modifications.  If mass_format is
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
                                        int precision);
/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the mass change of
 * any modifications in brackets following the modified residue.  If
 * If is mass_format MOD_MASSES_SEPARATE, all masses are listed in a
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
 MASS_FORMAT_T mass_format); // which mass value to print

/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the symbol for all
 * applied modifications.
 *
 * \returns A newly allocated array of characters, a text
 * representation of the modified sequence.
 */
char* modified_aa_string_to_string_with_symbols(
 MODIFIED_AA_T* aa_string, // modified aa's to translate
 int length); // length of aa_string

/**
 * \brief Takes an array of MODIFIED_AA_T's and returns an array of
 * char's, one for each amino acid in the sequence.  It DOES NOT
 * include any modification symbols.  Use with caution.
 * \returns A newly allocated char* with ONLY amino acids, all
 * modifications are removed.
 */
char* modified_aa_to_unmodified_string(MODIFIED_AA_T* aa_string, int length);

/**
 * \brief Allocates an array of MODIFIED_AA_T's and populates it with
 * the MODIFIED_AA_T value that corresponds to each sequence char
 * value and trailing modification symbols or masses.  Returns the new
 * sequence via the mod_sequence argument.
 *
 * \returns The length of the mod_sequence array.
 */
int convert_to_mod_aa_seq(const char* sequence, MODIFIED_AA_T** mod_sequence,
                          MASS_FORMAT_T mass_format = MOD_MASS_ONLY);

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(MODIFIED_AA_T* source, int length);

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(
  MODIFIED_AA_T* source ///< Sequence to copy
  );

/**
 * \returns whether the two modified sequences are equal or not
 */
bool equal_seq(
  const MODIFIED_AA_T* seq1, ///< Sequence 1
  const MODIFIED_AA_T* seq2 ///< Sequence 2
  );

/**
 * \brief Remove any characters not A-Z from a peptide sequence.
 * \returns A newly allocated string with the given sequence less any
 * modififcation symbols or masses.
 */
char* unmodify_sequence(const char* modified_seqeunce);

/**
 * \brief Remove any characters not A-Z from a peptide sequence.
 */
void unmodify_sequence_in_place(char* modified_seqeunce);

/**
 * \brief Determine if an array of MODIFIED_AA_T is a palindrome.  
 * Used by reverse_sequence to avoid returning a reversed sequence
 * that is the same as the target.  Ignores the first and last
 * residues. 
 * \returns TRUE if the reversed sequence would be the same as the
 * forward, otherwise FALSE.
 */
bool modified_aa_seq_is_palindrome(MODIFIED_AA_T* seq, int length);

/**
 * \brief Frees memory for an array of MODIFIED_AA_Ts.  Assumes is
 * terminated with the MOD_SEQ_NULL value
 */
void free_mod_aa_seq(MODIFIED_AA_T* seq);

/**
 * \brief Gives the size of the aa_mod struct.  For serialization
 */
int get_aa_mod_sizeof();

/**
 * The new definition of a PEPTIDE_T object.
 * 
 */
/*
struct peptide{
  unsigned char length; ///< The length of the peptide
  FLOAT_T peptide_mass;   ///< The peptide's mass with any modifications
  PEPTIDE_SRC_T* peptide_src; ///< a linklist of peptide_src
  bool is_modified;   ///< if true sequence != NULL
  MODIFIED_AA_T* sequence; ///< sequence with modifications
};
*/

/**
 * \brief checks to see if an amino acid is modified by a given mod
 * \returns TRUE if aa is modified by mod
 */
bool is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod);

/**
 * \brief Determine if this modified amino acid can be modified by
 * this modification.
 *
 * Checks the mod list of modifiable residues to see if aa is in the
 * list.  Also checks to see if aa has already been modified by this
 * mod.  
 * \returns TRUE if it can be modified, else FALSE
 */
bool is_aa_modifiable(MODIFIED_AA_T aa, AA_MOD_T* mod);

/**
 * \brief Adds a modification to a MODIFIED_AA_T.
 *
 * Assumes that the aa is modifiable, no explicit check.  If the aa is
 * already modified for the mod, no change to aa.
 */
void modify_aa(MODIFIED_AA_T* aa, const AA_MOD_T* mod);

/**
 * \brief Return the AA_MOD_T associated with the given symbol.  If
 * the symbol does not represent a modification, returns null.
 * Requires that parameters have been initialized.
 */
const AA_MOD_T* get_aa_mod_from_symbol(const char symbol);

/**
 * \brief Return the delta mass associated with the given modification
 * symbol.  If the symbol does not represent a modification, returns
 * 0. Requires that parameters have been initialized.
 */
FLOAT_T get_mod_mass_from_symbol(const char symbol);

/**
 * \brief Return the AA_MOD_T associated with the given mass shift.
 * The mass may either be from a single AA_MOD_T as given by the user
 * or from any combination of AA_MOD_T's.  If no AA_MOD_T(s) can be
 * found for the mass, returns null.
 * Requires that parameters have been initialized.
 */
const AA_MOD_T* get_aa_mod_from_mass(FLOAT_T mass);

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
bool compare_mods(AA_MOD_T** psm_file_mod_list, int num_mods);

/**
 * \brief Compare two mods to see if they are the same, i.e. same mass
 * change, unique identifier, position
 */
bool compare_two_mods(AA_MOD_T* mod1, AA_MOD_T* mod2);

/**
 * print all fields in mod.  For debugging
 */
void print_a_mod(AA_MOD_T* mod);

/* Setters and Getters */

/**
 * \brief Set the mass change caused by this modification.
 * \returns void
 */
void aa_mod_set_mass_change(AA_MOD_T* mod, double mass_change);
/**
 * \brief Get the mass change caused by this modification.
 * \returns The mass change caused by this modification.
 */
double aa_mod_get_mass_change(const AA_MOD_T* mod);

/**
 * \brief Access to the aa_list of the AA_MOD_T struct.  This pointer
 * can be used to get or set the list of residues on which this mod
 * can be placed.
 * \returns A pointer to the list of amino acids on which this mod can
 * be placed.
 */
bool* aa_mod_get_aa_list(AA_MOD_T* mod);

/**
 * \brief Set the maximum number of times this modification can be
 * placed on one peptide.
 * \returns void
 */
void aa_mod_set_max_per_peptide(AA_MOD_T* mod, int max);
/**
 * \brief Get the maximum number of times this modification can be
 * placed on one peptide.  
 * \returns The max times per peptide this mod can be placed.
 */
int aa_mod_get_max_per_peptide(AA_MOD_T* mod);

/**
 * \brief Set the maximum distance from the protein terminus that the
 * modification can be placed.  Which terminus (C or N) is determined
 * by the position type.  To indicate no position restriction, set to
 * MAX_PROTEIN_SEQ_LENGTH. 
 * \returns void
 */
void aa_mod_set_max_distance(AA_MOD_T* mod, int distance);
/**
 * \brief Get the maximum distance from the protein end that the
 * modification can be placed.  Will be MAX_PROTEIN_SEQ_LENGTH if
 * position type is ANY_POSITION.
 * \returns Maximum distance from protein terminus at which mod can be
 * placed. 
 */
int aa_mod_get_max_distance(AA_MOD_T* mod);

/**
 * \brief Set the position type of an aa_mod.
 * \returns void
 */
void aa_mod_set_position(AA_MOD_T* mod, MOD_POSITION_T position);

/**
 * \brief Where in the peptide can the modification be placed.
 * \returns ANY_POSITION for standard mods; C_TERM or N_TERM for those
 * that can only be placed at the ends of the peptide.
 */
MOD_POSITION_T aa_mod_get_position(AA_MOD_T* mod);

/**
 * \brief Sets whether the modification can prevent cleavage.
 * \returns void
 */
void aa_mod_set_prevents_cleavage(AA_MOD_T* mod, bool prevents_cleavage);

/**
 * \brief gets whether the modification can prevent cleavage
 * \returns TRUE or FALSE
 */
bool aa_mod_get_prevents_cleavage(AA_MOD_T* mod);

/**
 * \brief Sets whether the modifications can prevent cross-linking.
 * \returns void
 */
void aa_mod_set_prevents_xlink(AA_MOD_T* mod, bool prevents_xlink);

/**
 * \brief gets whether the modification can prevent cross-linking.
 * \returns TRUE or FALSE
 */
bool aa_mod_get_prevents_xlink(AA_MOD_T* mod);

/**
 * \brief The character used to uniquely identify the mod in the sqt file.
 * \returns The character identifier.
 */
char aa_mod_get_symbol(const AA_MOD_T* mod);

/**
 * \brief The bitmask used to uniquely identify the mod.
 * \returns The short int bitmask used to identify the mod.
 */
int aa_mod_get_identifier(const AA_MOD_T* mod);

/**
 * \brief Generates a string representation of an aa_mod and returns a
 * pointer to that newly allocated string.
 */
char* aa_mod_to_string(AA_MOD_T* mod);

/**
 * \brief Create a string containing all of the amino acids that can
 * be modified by this aa_mod.  E.g. if S, T, and Y can be modified,
 * returns "STY".
 * \returns A newly allocated string.
 */
char* aa_mod_get_aa_list_string(AA_MOD_T* mod);

/**
 * Count the number of modified aas in the string.
 */
int count_modified_aas(MODIFIED_AA_T* seq);

/**
 * /returns the mass of the modified sequence
 */
FLOAT_T get_mod_aa_seq_mass(
  MODIFIED_AA_T* seq, ///< The modified sequence
  MASS_TYPE_T mass_type ///<  mono or average?
  );

#endif //MODIFICATION_FILE_H


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

