/**
 * \file PostProcessProtein.h
 * AUTHOR : Sean McIlwain
 * $Revision: 1.00 $
 * \brief Object for representing a protein that doesn't exist within the database
 */
#ifndef POSTPROCESSPROTEIN_H_
#define POSTPROCESSPROTEIN_H_

#include "Protein.h"

/**
 *  There are cases where the protein database will not be available for the post-process
 *  routines. This class acts as a way to get around that problem.  The main difficulty is
 *  that the Peptide object relies on the PeptideSrc object to print out the sequence.
 *  This PeptideSrc object contains a Protein object and the start index of the
 *  peptide's sequence within this protein.  Parsing the search results files post-process
 *  without the database, only has the peptide sequence and the protein ids that the peptide
 *  came from.  PostProcessProtein keeps track of the peptides associated with a protein by 
 *  maintaining a vector. findStart returns the index of the peptide in the vector.  If 
 *  the peptide doesn't exist in the vector, than findStart will add the peptide to the vector 
 *  and return the index.
 *  getSequence and getSequencePointer to return the correct peptide sequence
 *  by using the offset parameter.  This way we can still use the Peptide and PeptideSrc objects
 *  together to print out the peptide sequence and the proteins associated with that sequence
 */

/**
 *  TODO: support the flanking amino acids so that peptide->getFlankingAAs() works correctly.
 */
class PostProcessProtein : public Crux::Protein {

  std::vector<std::string> sequences_;  ///< sequences that we have seen so far.
  std::vector<std::string> prev_aas_; ///< previous amino acid to the sequence
  std::vector<std::string> next_aas_; ///< next amino acid to the sequence

 public:

  
  /**
   * Default constructor
   */
  PostProcessProtein();

  /**
   * Default destructor
   */
  virtual ~PostProcessProtein();

  /**
   * \returns the index of the sequence in the sequences_ vector and
   * returns the index.  If sequence is not found, then add to the
   * list of sequences and return the last index
   */
  virtual int findStart(
    std::string sequence,  ///< the sequence to find
    std::string prev_aa,   ///< the previous aa of the sequence
    std::string next_aa    ///< the next aa of the sequence
  );

  /**
   * \returns the ith sequence
   */
  virtual char* getSequence(
    int offset=0 ///< The offset (or sequence index) for the sequence
  );

  /**
   * \returns the ith sequence pointer
   */
  virtual char* getSequencePointer(
    int offset=0 ///< The offset (or sequence index) for the sequence
  );

  char getNTermFlankingAA(
    int offset=0 ///< The offset (or sequence index) for the AA
  );

  char getCTermFlankingAA(
    int offset=0 ///< The offset (or sequence index) for the AA
  );

  /**
   * \returns true indicating that this is a PostProcessProtein object
   */
  virtual bool isPostProcess();

  /**
   * ends the program since we can't get the protein length from this
   * post process protein object
   */
  virtual unsigned int getLength();

  virtual Database* getDatabase();
 
};


#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
