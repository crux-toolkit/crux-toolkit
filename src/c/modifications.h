/**
 * \file modifications.h
 * \brief Datatypes and methods for peptide modifications
 *
 * Two data structures define modifications.  The AA_MOD_T is the most
 * basic type.  It is the information provided by the user: mass
 * change caused by this mod, amino acids which may be modified in
 * this way, and the maximum number of this type of modification which
 * may occur on one peptide.  A collection of AA_MODs that may occur
 * on some peptide are represented as a PEPTIDE_MOD_T.  This stores
 * a list of AA_MODS and the net mass change experienced by the
 * peptide.  AA_MODs are instantiated once as given by the user.  All
 * possible PEPTIDE_MODs are calcualted once and reused for each
 * spectrum search.  One PEPTIDE_MOD corresponds to one mass window
 * that must be searched.
 * 
 */

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
  MOD_POSITION_T position; ///< where the mod can occur, pep or prot position
  int max_per_peptide; ///< the maximum number of mods per peptide
  char symbol;         ///< the character to represent the mod in sqt files
  int identifier;      ///< the offset/bitmask assigned to this mod for unique
                       //identification, used with MODIFIED_AA
};
typedef struct _aa_mod AA_MOD_T;

/**
 * \enum _mod_position (typedefed as MOD_POSITION_T)
 * \brief An indication of where an AA_MOD may occur within a peptide.
 * Default is ANY_POSITION.
 */
enum _mod_position{ 
  ANY_POSITION,   ///< occur anywhere in any peptide
  PEPTIDE_C_TERM, ///< only c-terminus of peptide, seq[0]
  PEPTIDE_N_TERM, ///< only n-terminus of peptide, seq[len]
  PROTEIN_C_TERM, ///< only the c-term of the first peptide of a protein
  PROTEIN_N_TERM}; ///< only the n-term of the last peptide of a protein
/**
 * \typedef MOD_POSITION_T
 * \brief The typedef of the indicator for where an amino acid
 * modification can occur within a peptide and/or protein.
 */
typedef enum _mod_position MOD_POSITION_T;

// Can you think of a better name for this?
typedef short MODIFIED_AA_T; ///< letters in the expanded peptide
                             ///alphabet, bits for mod1 mod2...aa
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
 * \struct _peptide_mod
 *  A collection of aa modifications that can occur on a single
 *  peptide.  Includes the list of aa_mods and the net mass
 *  change of all mods.  Each aa_mod is included in the list as many
 *  times at it appears in the collection. There is no information about
 *  which residues in a peptide are modified.
 */
struct _peptide_mod{
  double mass_change;      ///< the net mass change for the peptide
  int num_uniq_mods;       ///< the number of items in the list_of_mods
  AA_MOD_T* list_of_mods;  ///< the list of aa_mods in this peptide
};
typedef struct _peptide_mod PEPTIDE_MOD_T*;


/**
 * \brief Generate a list of all PEPTIDE_MODs that can be considered
 * given the list of AA_MODs provided by the user.
 *
 * This only needs to be called once for a search.  The list of
 * PEPTIDE_MODs can be reused for each spectrum.
 *
 * \returns The number of peptide_mods returned in the
 * peptide_mod_list argument.
 */
int generate_peptide_mod_list(
 AA_MOD_T* aa_mod_list,           ///< list of aa_mods to permute
 int num_aa_mods,                 ///< number of itmes in above list
 PEPTIDE_MOD_T** peptide_mod_list ///< return here the list of peptide_mods
);

/**
 * \brief Check a peptide sequence to see if the aa_mods in
 * peptide_mod can be applied. 
 *
 * Assumes that an amino acid can be modified by more than one aa_mod,
 * but not more than once by a single aa_mod as defined in modifiable().
 * \returns TRUE if the sequence can be modified, else FALSE
 */
bool is_peptide_modifiable( PEPTIDE_T* peptide,
                            PEPTIDE_MOD_T* peptide_mod);


/**
 * \brief Take a peptide and a peptide_mod and return a list of
 * modified peptides.
 *
 * The peptide_mod should be guaranteed to be able to be applied to
 * the peptide at least once.  A single amino acid can be modified
 * multiple times by different aa_mods, but only once for a given
 * aa_mod as defined in modifiable().  
 * 
 * \returns The number of modified peptides in the array pointed to by
 * modified_peptides. 
 */
int modify_peptide(PEPTIDE_T* peptide,
                   PEPTIDE_MOD_T* peptide_mod,
                   PEPTIDE_T** modified_peptides);

/**
 * The new definition of a PEPTIDE_T object.
 * 
 */
struct peptide{
  unsigned char length; ///< The length of the peptide
  float peptide_mass;   ///< The peptide's mass with any modifications
  PEPTIDE_SRC_T* peptide_src; ///< a linklist of peptide_src
  BOOLEAN_T is_modified;   ///< if true sequence != NULL
  MODIFIED_AA_T* sequence; ///< sequence with modifications
};

/**
 * \brief checks to see if an amino acid is modified by a given mod
 * \returns TRUE if aa is modified by mod
 */
BOOLEAN_T is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod);


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
int aa_modified_by(MODIFIED_AA_T aa, 
                   AA_MOD_T* possible_mods, 
                   AA_MOD_T** mod_list);

char modified_aa_to_char(MODIFIED_AA_T aa);
MODIFIED_AA_T char_aa_to_modified(char aa);

// in parameter.c
//global
AA_MOD_T* list_of_mods;
int num_mods;
{MAX_AA_MODS = 11};// instead of #define; forces typechecking and obeys scope
AA_MOD_T position_mods[NUM_POSITIONS]; // if one of each 
// OR
AA_MOD_T* position_mods;
int num_position_mods;

/**
 * \brief Read the paramter file and populate the static parameter
 * list of AA_MODS.
 *
 * Also updates the array of amino_masses.  Dies with an error if the
 * number of mods in the parameter file is greater than MAX_AA_MODS.
 * \returns The number of mods requested.
 */
int read_mods_from_file(char* param_file, AA_MODS_T**);

/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.
 * \returns The number of mods pointed to by mods
 */
int get_aa_mods(AA_MOD_T*** mods);

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
BOOLEAN_T compare_mods(AA_MOD_T* psm_file_mod_list, int num_mods);

/**
 * \brief Compare two mods to see if they are the same, i.e. same mass
 * change, unique identifier, position
 */
BOOLEAN_T compare_two_mods(AA_MOD_T* mod1, AA_MOD_T* mod2);

// in mass.c
/**
 * \brief Extends the amino_masses table to include all possible
 * modifications.  
 *
 * Gets list of mods from parameter.c.  Should this fill in values for
 * both average and monoisotopic masses? 
 */
void extend_amino_masses(void);





