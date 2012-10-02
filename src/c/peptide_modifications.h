/**
 * \file peptide_modifications.h
 * \brief Datatypes and methods for peptide modifications
 *
 * Two data structures define modifications.  The AA_MOD_T is the most
 * basic type.  It is the information provided by the user: mass
 * change caused by this mod, amino acids which may be modified in
 * this way, and the maximum number of this type of modification which
 * may occur on one peptide.  AA_MODs are defined in modifications.h.
 * A collection of AA_MODs that may occur
 * on some peptide are represented as a PEPTIDE_MOD_T.  This stores
 * a list of AA_MODS and the net mass change experienced by the
 * peptide.  PEPTIDE_MODs are defined here.
 * AA_MODs are instantiated once after parsing the parameter file.  All
 * possible PEPTIDE_MODs are calcualted once and reused for each
 * spectrum search.  One PEPTIDE_MOD corresponds to one mass window
 * that must be searched.
 * 
 * $Revision: 1.5 $
 */
#ifndef PEPTIDE_MODIFICATIONS_H
#define PEPTIDE_MODIFICATIONS_H

#include <assert.h>
#include "utils.h"
#include "linked_list.h"
#include "objects.h"
#include "modifications.h"

/**
 * \brief Allocate a PEPTIDE_MOD and set all fields to default values
 * (i.e. no modifications).
 * \returns A heap allocated PEPTIDE_MOD with default values.
 */
PEPTIDE_MOD_T* new_peptide_mod();

/**
 * \brief Allocate a new peptide mod and copy contents of given mod
 * into it.
 * \returns A pointer to a new peptide mod which is a copy of the
 * given one.
 */
PEPTIDE_MOD_T* copy_peptide_mod(PEPTIDE_MOD_T* original);

/**
 * \brief Free the memory for a PEPTIDE_MOD including the aa_list.
 */
void free_peptide_mod(PEPTIDE_MOD_T* mod);

/**
 * \brief Generate a list of all PEPTIDE_MODs that can be considered
 * given the list of AA_MODs provided by the user and found in parameter.c.
 *
 * This only needs to be called once for a search.  The list of
 * PEPTIDE_MODs can be reused for each spectrum.
 *
 * \returns The number of peptide_mods returned in the
 * peptide_mod_list argument.
 */
int generate_peptide_mod_list(
 PEPTIDE_MOD_T*** peptide_mod_list ///< return here the list of peptide_mods
);

/**
 * \brief Check a peptide sequence to see if the aa_mods in
 * peptide_mod can be applied. 
 *
 * Assumes that an amino acid can be modified by more than one aa_mod,
 * but not more than once by a single aa_mod as defined in modifiable().
 * \returns TRUE if the sequence can be modified, else FALSE
 */
bool is_peptide_modifiable( Crux::Peptide* peptide,
                            PEPTIDE_MOD_T* peptide_mod);


/**
 * \brief Take a peptide and a peptide_mod and return via parameters a
 * list of modified peptides.
 *
 * The peptide_mod should be guaranteed to be able to be applied to
 * the peptide at least once.  A single amino acid can be modified
 * multiple times by different aa_mods, but only once for a given
 * aa_mod as defined in modifiable().  
 * 
 * \returns The number of modified peptides in the array pointed to by
 * modified_peptides. 
 */
int modify_peptide(Crux::Peptide* peptide,
                   PEPTIDE_MOD_T* peptide_mod,
                   LINKED_LIST_T* modified_peptides,
                   int max_aas_modified
);

/* Setters and Getters */
/**
 * \brief Add a new aa_mod to the peptide mod.  Updates mass_change,
 * num_mods and list_of_aa_mods.  Does not enforce the copy number of
 * an aa_mod to be less than max_per_peptide.
 * \returns void
 */
void peptide_mod_add_aa_mod(PEPTIDE_MOD_T* pep_mod,
                            //AA_MOD_T* aa_mod,
                            int aa_mod_idx,
                            int copies );

/**
 * \brief Get the value of the net mass change for this peptide_mod.
 * \returns The mass change for the peptide mod.
 */
double peptide_mod_get_mass_change(PEPTIDE_MOD_T* mod);

/**
 * \brief Get the number of aa_mods in this peptide_mod.
 * \returns The number of aa_mods in this peptide_mod.
 */
int peptide_mod_get_num_aa_mods(PEPTIDE_MOD_T* mod);

/**
 * \brief Get a pointer to the list of aa_mods in this peptide_mod.
 * The number of elements in the list is given by
 * peptide_mod_get_num_aa_mods. A unique aa_mod may be listed more
 * than once.  There is no particular order to the aa_mods in the
 * list.
 * \returns A pointer to a list of AA_MOD_T pointers.
 */
LINKED_LIST_T* peptide_mod_get_aa_mod_list(PEPTIDE_MOD_T* mod);

/**
 * \brief Compares the number of aa mods in two peptide mods for
 * sorting.
 * \returns Negative int, 0, or positive int if the number of aa_mods
 * in pmod 1 is less than, equal to or greater than the number of
 * aa_mods in pmod2, respectively.
 */
int compare_peptide_mod_num_aa_mods(const void* pmod1, 
                                    const void* pmod2);

/**
 * print all fields in peptide mod. For debugging
 */
void print_p_mod(PEPTIDE_MOD_T* mod);


#endif // PEPTIDE_MODIFICATIONS_H
















