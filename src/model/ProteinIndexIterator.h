#ifndef PROTEIN_INDEX_ITERATOR_H
#define PROTEIN_INDEX_ITERATOR_H
/**
 * \file ProteinIndexIterator.h
 * $Revision: 1.3 $
 * \brief Object for iterating a protein index file
 *****************************************************************************/

#include "Protein.h"

/**
 * protein index iterator
 * the protein index iterator parses the protein index file
 * for each protein index one by one and returns a new protein_index object
 */

class ProteinIndexIterator {
 protected:
  FILE* file_;  ///< The file handler of the fasta file
  Crux::Protein* next_protein_; ///< the next protein index to return
  bool has_next_; ///< is there a new protein to return?

  /**
   *
   *\returns true if successfully sets the ProteinIndexIterator, else false
   */
  bool setup();

 public:
    
  /**
   * 
   *\returns a new heap allocated protein index iterator
   */
  ProteinIndexIterator(
    const char* fasta_file ///< input fasta file -in
  );

  /**
   * Frees the allocated protein index iterator
   */
  virtual ~ProteinIndexIterator();

  /**
   *
   *\returns true if there is another protein index to return, else false
   */
  bool hasNext();

  /**
   *
   *\return the next protein index in the protein index file
   */
  Crux::Protein* next();

};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
