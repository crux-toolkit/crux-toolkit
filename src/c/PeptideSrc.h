/**
 * \file PeptideSrc.h
 * $Revision: 1.12 $
 * \brief Object for mapping a peptide to it's parent protein.
 */
#ifndef PEPTIDE_SRC_H
#define PEPTIDE_SRC_H
#include <stdio.h>

#include "utils.h"
#include "mass.h"

#include "objects.h"
#include "carp.h"
#include "PeptideConstraint.h"


class PeptideSrc {

 protected:
  DIGEST_T digestion_; ///< how specific the ends are relative to the enzyme
  Protein* parent_protein_; ///< the parent of this preptide
  int start_idx_; ///< start index of the peptide in the protein sequence, first residue is 1 
  PeptideSrc* next_association_; ///< a linklist of peptide_src     

  /**
   * Initializes an (empty) PeptideSrc object
   */
  void init();

 public:

  /**
   * \returns An (empty) peptide_src object.
   */
  PeptideSrc();

  /**
   * \returns a PROTEIN_PEPTIDE_ASSOCIATION object, populated with user
   * specified parameters 
   */
  PeptideSrc(
    DIGEST_T digest,
    Protein* parent_protein, ///< the parent of this preptide -in
    int start_idx ///< peptide start index in protein sequence, first is 1 -in
    );

  /**
   *\returns an array of PROTEIN_PEPTIDE_SRC object
   * only used in index.c, when the peptide src count for  peptide is known
   */
  static PeptideSrc* newArray(
    int size ///< the size of the peptide_src array -in
    );

  /**
   * \brief Fill in the values from the original array into the new array.
   * Assumes that the new array has been allocated by newArray().
   */
  static void copyArray(
    PeptideSrc* original_array, 
    PeptideSrc* new_array, 
    int array_size);

  /**
   *\returns a linklist of PROTEIN_PEPTIDE_SRC object
   * only used in index.c, when the peptide src count for peptide is known
   */
  static PeptideSrc* newLinklist(
    int size ///< the size of the peptide_src array -in
    );

  /**
   *\returns the PROTEIN_PEPTIDE_SRC object in the array with the index
   * index starts at 0.
   * only used in index.c, when the peptide src count for  peptide is known
   */
  static void setArray(
    PeptideSrc* src_array , ///< the working peptide src_arry -out
    int array_idx, ///< array index of the peptide_src to set
    DIGEST_T digest,
    Protein* parent_protein, ///< the parent of this preptide -in
    int start_idx ///< start index of the peptide in the protein sequence -in
    );

  /**
   * Frees the entire allocated peptide_src linklist object
   * Assumes that peptide src is Link list implementation
   */
  static void free(
    PeptideSrc* peptide_src ///< object to free -in 
    );

  /**
   * Frees the an individual allocated peptide_src object
   * assumes that new_association pointer is NULL or some other pointer exist for the rest of the linklist 
   */
  virtual ~PeptideSrc();


  /**
   * Prints a peptide object to file.
   */
  void print(
    FILE* file  ///< the output stream -out
    );

  /**
   * \brief Read in the peptide_src objects from the given file and
   * assosiated them with the given peptide.  
   * Proteins for the pepitde_src are found in the given database.  If
   * database is NULL, does not set proteins.  (This option is used for
   * sorting index files while creating index.)  Either array or 
   * linked list implementation of multiple peptide_src is used based on
   * the value of use_array.
   *
   * \returns true if peptide_src's were successfully parsed, else
   * returns false.
   */
  static bool parseTabDelimited(
    Peptide* peptide,   ///< assign peptide_src(s) to this peptide
    MatchFileReader& file,           ///< file to read from
    Database* database, ///< database containing proteins
    bool use_array); ///< use array implementation vs. linked list

  /**
   * \brief Read in the peptide_src objects from the given file and
   * assosiated them with the given peptide.  
   *
   * Proteins for the pepitde_src are found in the given database.  If
   * database is NULL, does not set proteins.  (This option used for
   * sorting index files while creating index.) Either array or linked
   * list implementation of multiple peptide_src is used based on the
   * value of use_array. 
   *
   * \returns true if peptide_src's were successfully parsed, else
   * returns false and sets peptide's peptide_src member variable to
   * NULL. 
   */
  static bool parse(
    Peptide* peptide,   ///< assign peptide_src(s) to this peptide
    FILE* file,           ///< file to read from
    Database* database, ///< database containing proteins
    bool use_array);///< use array implementation vs. linked list

  /**
   * Copies the entire linklist of peptide_src object src to dest.
   * dest must be a heap allocated peptide_src
   */
  static void copy(
    PeptideSrc* src, ///< source peptide_src -in
    PeptideSrc* dest ///< destination peptide_src -out
    );

  /**
   * sets the level of digestion
   */
  void setDigest(
    DIGEST_T digest ///< the type of the peptide -in
    );

  /**
   * \returns the level of digestion
   */
  DIGEST_T getDigest();

  /**
   * sets the parent protein
   */
  void setParentProtein(
    Protein* parent_protein ///< the parent of this preptide -in  
    );

  /**
   * \returns a pointer to the parent protein
   */
  Protein* getParentProtein();

  /**
   * sets the start index of the peptide in the protein sequence
   */
  void setStartIdx(
    int start_idx ///< start index of the peptide in the protein sequence -in
    );

  /**
   * \returns the start index of the peptide in the protein sequence
   */
  int getStartIdx();

  /**
   * sets the next peptide_src on the link list
   * assumes that the src_association's next_association field is NULL
   */
  void setNextAssociation(
    PeptideSrc* new_association ///< the new peptide_src to add -in   
    );

  /**
   * \returns the next peptide_src on the link list
   */
  PeptideSrc* getNextAssociation();

  /**
   * \returns a pointer to the start of the peptide with in it's parent protein sequence
   */
  char* getSequencePointer();

  /**
   *\returns the peptide_src strct size, value of sizeof function
   */
  static int getSizeOf(void);

  /**
   * serialize peptide src in binary
   */
  void serialize(
  FILE* file  ///< output file -in   
  );

  /**
   * Return the number of bytes taken up by one peptide_src when
   * serialized to file.  Used for skipping past peptide_src in an index
   * file. 
   */
  static int sizeOfSerialized();
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
