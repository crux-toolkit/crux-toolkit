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
#include "Protein.h"
#include "objects.h"
#include "carp.h"
#include "PeptideConstraint.h"
#include "Database.h"
#include "modifications.h"
#include "peptide_modifications.h"

//these may be elsewhere
static const int MAX_PEPTIDE_LENGTH = 255;

#include <string>

/**
 * \class peptide
 * \brief A subsequence of a protein.
 */

namespace Crux {

class Peptide {

 protected:

  /**
   * static global variable
   * determines if the peptide src are created by link lists or array
   * if true, peptides are implented with link list peptide src, else array
   */

  /* Private data types */
  unsigned char length_; ///< The length of the peptide
  FLOAT_T peptide_mass_;   ///< The peptide's mass.
  std:: vector<PeptideSrc*> peptide_srcs_; ///< a vector of peptide_srcs_

  MODIFIED_AA_T* modified_seq_; ///< peptide sequence with modifications
  MODIFIED_AA_T* decoy_modified_seq_; ///< randomized peptide sequence

  /**
   * Initializes an (empty) peptide object
   */
  void init();

 public:

  /*  Allocators/deallocators  */
  
  /**
   * \returns An (empty) peptide object.
   */
  Peptide();

  /**
   * \returns A new peptide object, populated with the user specified
   * parameters.
   */
  Peptide(
    unsigned char length,     ///< The length of the peptide -in
    FLOAT_T peptide_mass,       ///< The neutral mass of the peptide -in
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
   *\returns the protein struct size, value of sizeof function
   */
  int getSizeOf(void);
 
  /**
   * Merge to identical peptides, copy all peptide_src into one of the peptide
   * peptide_dest, peptide_bye must have at least one peptide src
   * frees the peptide_bye, once the peptide_src are re-linked to the peptide_dest
   * Assumes that both peptides use linklist implemenation for peptide_src
   * \returns true if merge is successful else false
   */
  static bool mergePeptides(
    Peptide* peptide_dest,
    Peptide* peptide_bye
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
  
  /*  Get-set:  mass */


  /**
   * sets the peptide mass
   */
  void setPeptideMass(
    FLOAT_T peptide_mass  ///< the mass of the peptide - in
    );

  /**
   * \returns the peptide mass
   */
  /*inline*/ FLOAT_T getPeptideMass();

  /** 
   * \returns the mass of the peptide if it had charge "charge"
   */
  FLOAT_T getChargedMass(
    int charge ///< charge of peptide -in
    );

  /** 
   * \returns the m/z of the peptide if it had charge "charge"
   */
  FLOAT_T getMz(
    int charge ///< the charge of peptide -in
    );

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
  PeptideSrc* getPeptideSrc();

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

  /*  Get-set:  sequence */

  /**
   * sets the sequence length of the peptide
   */
  void setLength(
    unsigned char length  ///< the length of sequence -in
    );

  /**
   *\returns the sequence length of the peptide
   */
  unsigned char getLength();

  /**
   * \brief Get the sequence of a peptide.
   * Goes to the first peptide_src to gain sequence, thus must have at
   * least one peptide src 
   * \returns A newly allocated copy of the sequence.
   */
  char* getSequence();

  /**
   * \brief Get a string representation of the target (unshuffled)
   * peptide sequence with no added modification symbols.
   * For target peptides, returns the same as get_peptide_sequence.
   * \returns The newly-allocated sequence of peptide
   */
  char* getUnshuffledSequence();

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
  char* getSequenceSqt();

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
  char* getSequenceFromPeptideSrcSqt(
    PeptideSrc* peptide_src ///< peptide_src -in 
   );

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

  /**
   * \brief Get the modified aa sequence
   *
   * If the peptide has no modifications, create a sequence of
   * MODIFIED_AA_T's in which none of them are actually modified.
   * \returns A newly allocated copy of the sequence of MODIFIED_AA_Ts.
   */
  //MODIFIED_AA_T* get_peptide_modified_sequence( // why is this not working??!!
  unsigned short* getModifiedAASequence();
  
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
   * \returns A newly allocated string of the peptide sequence including
   * any modifications.
   */
  char* getModifiedSequenceWithSymbols();

  /**
   * \brief Get the modified aa sequence in string form.
   *
   * If the peptide has no modifications, returns same string as
   * get_peptide_sequence.  If modified, adds in brackets the masses of
   * all modifications.  If merge_masses is true, prints the sum of all
   * modifications for a residue.  If false, prints all masses in a
   * comma separated list.
   * \returns A newly allocated string of the peptide sequence including
   * any modifications.
   */
  char* getModifiedSequenceWithMasses(
    MASS_FORMAT_T merge_masses ///< do we want to merge masses?
    );

  /**
   * \brief Get the target sequence of the peptide encoded as char*
   * including modification symbols (e.g. *,#).
   *
   * If the peptide is not a decoy, returns the same sequence as
   * get_peptide_modified_sequence.  If the peptide has no
   * modifications, returns same string as get_peptide_sequence.  If
   * modified, adds the mod symbols to the string. 
   * \returns A newly allocated string of the peptide's unshuffled
   * (target) sequence including any modifications.
   */
  char* getUnshuffledModifiedSequence();

  /*  Getters requiring calculation */
  int countModifiedAAs();

  /**
   * \returns The mass of the given peptide.
   */
  static FLOAT_T calcSequenceMass(
    const char* peptide, ///< the query peptide -in
    MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
    );

  /**
   * \returns The mass of the given peptide.
   */
  FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
    );

  /**
   * \returns The hydrophobicity of the given peptide, as in Krokhin (2004).
   */
  FLOAT_T calcKrokhinHydrophobicity();

  /**
   * Examines the peptide sequence and counts how many tryptic missed
   * cleavage sites exist. 
   *\returns the number of missed cleavage sites in the peptide
   */
  int getMissedCleavageSites();

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
   * Creates a heap allocated hash_value for the peptide that should
   * uniquely identify the peptide
   *\returns the string of "<first src protein idx><start idx><length>"
   */
  char* getHashValue();

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

  /*  Comparisons for sorting  */
  
  /**
   * Compare peptide sequence
   * \returns true if peptide sequence is identical else false
   */
  static bool compareSequence(
    Peptide* peptide_one,
    Peptide* peptide_two
  );

  /**
   * Compare two peptide sequences.
   * \returns Zero (0) if the sequences are identical, -1 if the first
   * sequence is less than the first and 1 if the first sequence is
   * greater than teh first.
   */
  static int triCompareSequence(
    Peptide* peptide_one,  ///< the peptide sequence to compare  -out
    Peptide* peptide_two  ///< the peptide sequence to compare  -out
    );

  /**
   * Compare the sequence of two peptides and return true if the first
   * petpide sequence is less than (in a lexical sort) the second peptide.
   */
  static bool lessThan(
    Peptide* peptide_one,
    Peptide* peptide_two
    );

  /**
   * compares two peptides with the lexical sort type
   * for qsort
   * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
   */
  static int compareLexicalQSort(
    Peptide** peptide_one, ///< peptide to compare one -in
    Peptide** peptide_two ///< peptide to compare two -in
    );

  /**
   * compares two peptides with the mass sort type
   * if peptide mass is identical sort by lexicographical order
   * used for qsort function
   * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
   */
  static int compareMassQSort(
    Peptide** peptide_one, ///< peptide to compare one -in
    Peptide** peptide_two ///< peptide to compare two -in
    );

  /**
   * compares two peptides with the length sort type
   * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
   */
  static int compareLengthQSort(
    Peptide** peptide_one, ///< peptide to compare one -in
    Peptide** peptide_two ///< peptide to compare two -in
    );

  /**
   * Compare peptide mass
   * \returns 0 if peptide mass is identical else 1 if peptide_one is larger, -1 if peptide_two is larger
   */
  static int compareMass(
    Peptide* peptide_one,
    Peptide* peptide_two
    );

  /*  Printing / parsing       */

  /**
   * Prints a peptide object to file.
   * prints all peptide_src object it's associated 
   * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
   *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-trypticity> <\\t peptide-sequence> \n
   * prints in correct format for generate_peptide
   */
  void printInFormat(
    bool flag_out, ///< print peptide sequence? -in
    FILE* file  ///< the out put stream -out
    );

  /**
   * Prints a peptide object to file.
   * ONLY prints peptide_src that match the peptide_src
   * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
   *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
   * prints in correct format for generate_peptide
   */
  void printFilteredInFormat(
    bool flag_out, ///< print peptide sequence? -in
    FILE* file  ///< the out put stream -out
    );

  /**
   * Serialize a peptide to a FILE in binary
   * \returns true if serialization is successful, else false
   *
   * The peptide serialization format looks like this:
   *
   *<Peptide: peptide struct><int: number of peptide_src>[<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>]+
   * the bracket peptide src information repeats for the number of peptide src listed before the bracket
   * the protein index is the index of the parent protein in the database Database
   *
   */
  bool serialize(
    FILE* file,
    FILE* text_file
    );
 
  /**
   * \brief Read in a peptide from a tab-delimited file and return it.
   *
   * Parses the information for a peptide match from the search
   * file.  Allocates memory for the peptide and all of
   * its peptide_src's.  Requires a database so that the protein can be
   * set for each peptide_src.  Returns NULL if eof or if file format
   * appears incorrect.
   *
   * \returns A newly allocated peptide or NULL
   */
  static Peptide* parseTabDelimited(
    MatchFileReader& file, ///< the tab delimited peptide file -in
    Database* database,///< the database containing the peptides -in
    Database* decoy_database = NULL ///< optional database with decoy peptides
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
  static Peptide* parse(
    FILE* file, ///< the serialized peptide file -in
    Database* database ///< the database containing the peptides -in
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
   * \returns true if peptide was successfully parsed or false if it was
   * not. 
   */
  bool parseNoSrc(
    FILE* file,       ///< file pointing to a serialized peptide
    long int* pepitde_src_file_location);  // use to seek back to peptide_src

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
   *
   * \returns a pointer to the string. Caller is responsible for freeing memeory.
   * If peptide has no sources returns NULL.
   */
  char *getProteinIds();

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
