/**
 * \file peptide.h 
 * $Revision: 1.52 $
 * \brief Object for representing one peptide.
 */
#ifndef CRUX_PEPTIDE_H_
#define CRUX_PEPTIDE_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util/utils.h"
#include "util/crux-utils.h"
#include "util/mass.h"
#include "Modification.h"
#include "Protein.h"
#include "model/objects.h"
#include "io/carp.h"
#include "PeptideConstraint.h"
#include "Database.h"
#include "util/modifications.h"
#include "util/peptide_modifications.h"

//these may be elsewhere
static const int MAX_PEPTIDE_LENGTH = 255;

#include <string>

/**
 * \class peptide
 * \brief A subsequence of a protein.
 */

class ModificationDefinition;

namespace Crux {

class Modification;

class Peptide {

 protected:
  unsigned char length_; ///< The length of the peptide
  std::vector<PeptideSrc*> peptide_srcs_; ///< a vector of peptide_srcs_

  MODIFIED_AA_T* decoy_modified_seq_; ///< randomized peptide sequence
  std::string sequence_;
  std::vector<Modification> varMods_;

 public:
  /*  Allocators/deallocators  */
  
  /**
   * \returns An (empty) peptide object.
   */
  Peptide();
  Peptide(std::string sequence);
  Peptide(std::string sequence, std::vector<Modification> mods);

  /**
   * \returns A new peptide object, populated with the user specified
   * parameters.
   */
  Peptide(
    unsigned char length,     ///< The length of the peptide -in
    Crux::Protein* parent_protein, ///< The parent_protein of this peptide -in
    int start_idx ///< Start index of peptide in the protein sequence -in
    );

  /**
   * \brief Allocates a new peptide giving it the values of the source
   * peptide.
   * \returns A newly allocated peptide identical to the source.
   */
  Peptide(
    Peptide* src ///< source peptide -in
  );
                             
  /**
   * Merges two identical peptides by adding the peptide_src of the
   * second to the first.  The second peptide remains unchanged.
   * Does not comfirm identity of peptides.
   * \returns true if merge is successfull.
   */
  static bool mergePeptidesCopySrc(
    Peptide* peptide_dest,
    Peptide* peptide_giver
    );

  /**
   * Frees an allocated peptide object.
   * Depending on peptide_src implementation determines how to free srcs
   * This decision is made by global variable PEPTIDE_SRC_USE_LINK_LIST
   */
  ~Peptide();

  /*  Getters and Setters  */

  /*  Get-set:  source */
  
  /**
   * sets the peptide_src field in the peptide
   * must pass on a heap allocated peptide_src object
   * does not copy in the object, just the pointer to the object.
   */
  void setPeptideSrc(
    PeptideSrc* new_association ///< new peptide_src -in
    );

  /**
   * this method adds the new_association to the end of the existing peptide's 
   * linklist of peptide_srcs
   * must pass on a heap allocated peptide_src object
   * does not copy in the object, just the pointer to the object.
   */
  void addPeptideSrc(
    PeptideSrc* new_association ///< new peptide_src -in
    );

  /**
   * this method adds the peptide src array to an EMPTY peptide
   * only used in index.c, when the peptide src count for  peptide is known
   * Any existing peptide_src will lose it's reference
   */
  void addPeptideSrcArray(
    PeptideSrc* peptide_src_array ///< new peptide_src -in
    );

  /**
   * returns a pointer to the first PeptideSrc object of the peptide
   */
  PeptideSrc* getPeptideSrc() const;

  /**
   * returns a point to the peptide_protein_association field of the peptide
   */
  std::vector<PeptideSrc*>& getPeptideSrcVector();

  PeptideSrcIterator getPeptideSrcBegin();
  PeptideSrcIterator getPeptideSrcEnd();

  /**
   * get the peptide->first peptide_src->parent protein->database
   */
  Database* getFirstSrcDatabase();

  /**
   * \returns The number of peptide sources (i.e. proteins) the peptide has.
   */
  int getNumPeptideSrc();

  /**
   * returns a pointer to the peptide's first parent protein field of the peptide
   */
  Crux::Protein* getParentProtein();

  /**
   * sets the sequence length of the peptide
   */
  void setLength(
    unsigned char length  ///< the length of sequence -in
    );

  /**
   *\returns the sequence length of the peptide
   */
  unsigned char getLength() const;

  /**
   * \brief Get the sequence of a peptide.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   * \returns A newly allocated copy of the sequence.
   */
  char* getSequence() const;

  /**
   * \brief Get a string representation of the target (unshuffled)
   * peptide sequence with no added modification symbols.
   * For target peptides, returns the same as get_peptide_sequence.
   * \returns The newly-allocated sequence of peptide
   */
  std::string getUnshuffledSequence() const;

  /**
   * \returns a pointer to the start of peptide sequence with in it's protein parent sequence, 
   * thus does not have terminating signe until end of parent protein
   * goes to the first peptide_src to find the location of start, thus must have at least one peptide src
   * should not print, will result in printing the entire protein sequence
   */
  char* getSequencePointer();

  /**
   * \brief Formats the sequence of the peptide with each flanking AA.
   * 
   * Format is "X.peptide_sequence.X", where "X" is a flanking amino acid.
   * "X", is printed as "-" if there is no flanking sequence.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   * \returns A newly allocated string with the sqt-formated peptide sequence.
   */
  std::string getSequenceSqt();

  /**
   * \brief Return a char for the amino acid c-terminal to the peptide
   * in the peptide src at the given index.
   *
   * \returns A char (A-Z) or - if peptide is the first in the protein.
   */
  char getCTermFlankingAA();

  /**
   * \brief Return a char for the amino acid n-terminal to the peptide
   * in the peptide src at the given index.
   *
   * \returns A char (A-Z) or - if peptide is the last in the protein.
   */
  char getNTermFlankingAA();

  void addMod(const ModificationDefinition* mod, unsigned char index);
  void setMods(const std::vector<Modification>& mods);
  std::vector<Modification> getMods() const;
  std::vector<Modification> getVarMods() const;
  std::vector<Modification> getStaticMods() const;

  bool hasMonoLink() const;
  
  /**
   * \brief Add a modification to a peptide.
   *
   * Adds the modified sequence to the peptide and changes the peptide
   * mass based on the mass change in the peptide_mod.
   * \returns void
   */
  void setMod(
    MODIFIED_AA_T* mod_seq, ///< modified seq to add
    PEPTIDE_MOD_T* pep_mod  ///< mod that made the seq
  );

  std::string getModsString() const;

  bool isModified();

  bool isDecoy();

  std::string getDecoyType();

  /**
   * \brief Get the modified aa sequence
   *
   * If the peptide has no modifications, create a sequence of
   * MODIFIED_AA_T's in which none of them are actually modified.
   * \returns A newly allocated copy of the sequence of MODIFIED_AA_Ts.
   */
  //MODIFIED_AA_T* get_peptide_modified_sequence( // why is this not working??!!
  unsigned short* getModifiedAASequence();

  static std::string unmodifySequence(const std::string& seq);
  void setUnmodifiedSequence(const std::string& sequence);
  
  /**
   * sets the modified sequence for the peptide
   */
  void setModifiedAASequence(
    MODIFIED_AA_T* mod_seq,  ///< modified sequence to set
    bool decoy ///< is the peptide a decoy?
  );

  /**
   * \brief Get the modified aa sequence in string form.
   *
   * If the peptide has no modifications, returns same string as
   * get_peptide_sequence.  If modified, adds the mod symbols to the string.
   * \returns The peptide sequence including any modifications.
   */
  std::string getModifiedSequenceWithSymbols();

  /**
   * \brief Get the modified aa sequence in string form.
   *
   * If the peptide has no modifications, returns same string as
   * get_peptide_sequence.  If modified, adds in brackets the masses of
   * all modifications.  If merge_masses is true, prints the sum of all
   * modifications for a residue.  If false, prints all masses in a
   * comma separated list.
   * \returns The peptide sequence including any modifications.
   */
  std::string getModifiedSequenceWithMasses();

  void setDecoyModifiedSeq(MODIFIED_AA_T* decoy_modified_seq);

  /*  Getters requiring calculation */
  int countModifiedAAs();

  /**
   * \returns The mass of the given peptide.
   */
  static FLOAT_T calcSequenceMass(
    const std::string& peptide, ///< the query peptide -in
    MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
    );

  /**
   * \returns The mass of the given peptide.
   */
  FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
    ) const;

  /**
   * \returns the mass of the given peptide, with modifications
   */
  FLOAT_T calcModifiedMass(
    MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  ) const;

  FLOAT_T calcModifiedMass() const;

  /**
   * Examines the peptide sequence and counts how many tryptic missed
   * cleavage sites exist. 
   *\returns the number of missed cleavage sites in the peptide
   */
  int getMissedCleavageSites();

  int getMissedCleavageSites(
    std::set<int> skip //skip these amino acid indices.
  );

  /**
   * \brief Find the distance from the n-terminus of the source protein
   * to the n-terminus of the peptide.  
   * In the case of multiple source proteins, return the smallest
   * distance.
   * \returns The distance from the protein n-terminus.
   */
  int getNDistance();

  /**
   * \brief Find the distance from the c-terminus of the source protein
   * to the c-terminus of the peptide.
   * In the case of multiple source proteins, return the smallest
   * distance.
   * \returns The distance from the protein c-terminus.
   */
  int getCDistance();

  /**
   * Change the given target peptide into a decoy by randomizing its sequence.
   * Uses settings in parameter.c to decide between shuffling and
   * reversing the sequence.  Any modifications that exist will be
   * maintained on the same amino acids whose position will move.
   */
  void transformToDecoy();

  /**
   * \brief Return a randomly shuffled version of the given peptide's 
   * sequence as an array of char (A-Z).  Based on the peptide type,
   * will leave the end(s) unchanged to preserve the tryptic property. 
   * 
   *\returns A newly-allcoated char array of the shuffled sequence.
   */
  char* generateShuffledSequence();

  /**
   * \brief Return a reversed version of the given peptide's sequence as
   * an array of char (A-Z).  Leave the first and last residue
   * unchanged.  If the reversed sequence is identical to the target,
   * shuffle the sequence instead.
   *
   * \returns A newly-allocated char array of the reversed sequence.
   */
  char* generateReversedSequence();

  /**
   * \brief Return a randomly shuffled version of the given peptide's 
   * sequence as an array of MODIIFIED_AA_T.  Based on the peptide type,
   * will leave the end(s) unchanged to preserve the tryptic property.
   * 
   *\returns A newly-allcoated MODIFIED_AA_T array of the shuffled sequence.
   */
  MODIFIED_AA_T* generateShuffledModSequence();

  /**
   * \brief Return a reversed version of the given peptide's sequence as
   * an array of MODIFIED_AA_T.  Leave the first and last residue
   * unchanged.  If the reversed sequence is identical to the target,
   * shuffle the sequence instead.
   *
   * \returns A newly-allocated MODIFIED_AA_T array of the reversed sequence.
   */
  MODIFIED_AA_T* generateReversedModSequence();

  std::string getId() const;
  static std::string getId(const std::string& unmodifiedSeq, const std::vector<Modification>& mods);

  /*  Comparisons for sorting  */

  /**
   * Compare two peptide sequences.
   * \returns Zero (0) if the sequences are identical, -1 if the first
   * sequence is less than the first and 1 if the first sequence is
   * greater than the first.
   */
  static int triCompareSequence(
    Peptide* peptide_one,  ///< the peptide sequence to compare  -out
    Peptide* peptide_two  ///< the peptide sequence to compare  -out
    );

  /**
   * Compare the sequence of two peptides and return true if the first
   * peptide sequence is less than (in a lexical sort) the second peptide.
   */
  static bool lessThan(
    Peptide* peptide_one,
    Peptide* peptide_two
    );

  /**
   * \brief Builds a comma delimited string listing the 
   * protein id(peptide start index) for the sources of 
   * a peptide
   *
   * \returns a string of the protein sources for this peptide
   */
  std::string getProteinIdsLocations();

  /**
   * \brief Builds a comma delimited string listing the protein ids
   * for the sources of a peptide.
   */
  std::vector<std::string> getProteinIds();

  /**
   * \brief Builds a comma delimited string listing the flanking amino acids
   * for the sources of a peptide.
   *
   * \returns a pointer to the string. Caller is responsible for freeing memeory.
   * If peptide has no sources returns NULL.
   */
  char* getFlankingAAs();

  /**
   * Fills the given vectors with the names and descriptions of all
   * proteins containing this peptide.  Returned in the same order as
   * getFlankingAAs().  Clears any existing data in the vectors.
   * \returns The number of proteins.
   */
  int getProteinInfo(std::vector<std::string>& protein_ids,
                     std::vector<std::string>& protein_descriptions);

};  // class Peptide

};  // namespace Crux

/*  Iterators */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  Crux::Peptide* peptide ///< peptide sequence to iterate -in
  );

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  );

/**
 * The basic iterator functions.
 * \returns true if there are additional residues to iterate over, false if not.
 */
bool residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  );

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  );

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
