/**
 * \file modifications.c
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
 */

/**
 * \brief Generate a list of all PEPTIDE_MODs that can be considered
 * given the list of AA_MODs provided by the user.
 *
 * This only needs to be called once for a search.  The list of
 * PEPTIDE_MODs can be reused for each spectrum.
 *
 * \returns The number of peptide_mods listed in the
 * peptide_mod_list argument.
 */
int generate_peptide_mod_list(
 AA_MOD_T* aa_mod_list,           ///< list of aa_mods to permute
 int num_aa_mods,                 ///< number of itmes in above list
 PEPTIDE_MOD_T** peptide_mod_list ///< return here the list of peptide_mods
 ){
  /*
    get the number of items in final list

    init final_list
    init temp_list

    for each mod in aa_mod_list
      for i 1 to mod_max
        create a new peptide_mod with i mods
        add to the temp_list
        for each pep_mod in final_list
          copy the pep_mod and add i aa_mods
          add to the temp_list
        end of final_list
      end i==max
      add the temp_list to the final_list
    last mod
   */
{

/**
 * \brief Check a peptide sequence to see if the aa_mods in
 * peptide_mod can be applied. 
 *
 * Assumes that an amino acid can be modified by more than one aa_mod,
 * but not more than once by a sligle aa_mod.
 * \returns TRUE if the sequence can be modified, else FALSE
 */
bool is_peptide_modifiable( PEPTIDE_T* peptide,
                            PEPTIDE_MOD_T* peptide_mod){


  /* 

  for each aa_mod in peptide_mod
    initialize a counter to 0
  last aa_mod

  for each aa in peptide sequence
    for each aa_mod
      counter[aa_mod] += aa_mod.aa_list[aa]
    last aa_mod 
  last aa

  bool success = true
  for each aa_mod in peptide_mod
    success && (counter[aa_mod] >= peptide_mod.occurs[aa_mod])
  last aa_mod

  return success
  */
}

/**
 * \brief Take a peptide and a peptide_mod and return a list of
 * modified peptides.
 *
 * The peptide_mod should be guaranteed to be applicable to
 * the peptide at least once.  A single amino acid can be modified
 * multiple times by different aa_mods, but only once for a given
 * aa_mod (as in defined by modifiable()).  
 * 
 * \returns The number of modified peptides in the array pointed to by
 * modified_peptides. 
 */
int modify_peptide(
  PEPTIDE_T* peptide,             ///< the peptide to modify
  PEPTIDE_MOD_T* peptide_mod,     ///< the set of aa_mods to apply
  PEPTIDE_T** modified_peptides){ ///< the returned modified peptides

  /*
  MODIFIED_AA_T** modified_seqs = NULL
  initialize array to hold the peptide seq with the first aa_mod applied
  and set count to the size of the array
  (ie int count = apply_mods_to_seq(peptide_seq, aa_mod[1], &modified_seqs)

  for each aa_mod in peptide_mod
    count = apply_mods_to_list(&mod_pep_seqs, count, aa_mod)
  last aa_mod

  for each seq in mod_pep_seq
    copy peptide and give it the modified seq
    add to the return array
  last seq

  return count
   */
}


/**
 * take a list of modified_peptide_seqs, add to it the given
 * modification in all possible permutations
 */
int apply_mods_to_list(
  MODIFIED_AA_T*** mod_seqs, ///< a pointer to an array of arrays
  int count, ///< the number of arrays pointed to by mod_seqs*
  AA_MOD_T* mod_to_apply){  ///< the specific mod to apply

  /*
    init final list
    for each mod_seq
      init added_seqs
      count += apply_mods_to_seq(seq, mod, &added_seqs)
      add added_seqs to final list
    last mod_seq

    replace mod_seqs with final_list
    return count
   */
}

int apply_mods_to_seq(
  MODIFIED_AA_T* seq, ///< the seq to modify
  AA_MOD_T* mod,      ///< the mod to apply
  MODIFIED_AA_T*** return_list){ ///< the newly modified versions of
                                 ///seq

  /*
    init return_list
    for each aa in seq
      if modifiable(aa, mod)
      then 
         copy seq
         modify the aa in seq
         add to list
         increment count
      end if
    last aa

    return count
   */
  
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
bool modifiable(MODIFIED_AA_T* aa, AA_MOD_T* mod){
  // if the residue is in the list of residues that can be modified
  // AND if the residue is not already modified by the same mod
  //    return true
  // else 
  //    return false
}

BOOLEAN_T is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod){
  int id = mod->unique_id;
  if( (aa && id) != 0 ){
    return TRUE;
  } 
  return FALSE;
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
int aa_modified_by(MODIFIED_AA_T aa, 
                   int num_possible,
                   AA_MOD_T* possible_mods, 
                   AA_MOD_T*** mod_list){

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

}


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


