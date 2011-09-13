/**
 * \file ProteinIndex.h
 * $Revision: 1.3 $
 * \brief Object for creating a protein index
 *****************************************************************************/
#ifndef PROTEIN_INDEX_H
#define PROTEIN_INDEX_H

#include <stdio.h>
#include "utils.h"
#include "objects.h"
#include "Peptide.h"
#include "Protein.h"
#include "carp.h"
#include "PeptideConstraint.h"


class ProteinIndex {
 protected:
  unsigned long int offset_;  ///< The offset of the protein in the fasta file
  unsigned int protein_idx_;  ///< The idx of the protein in the fasta file
 public:


  /**
   * \returns An (empty) protein_index object.
   */
  void init();
  ProteinIndex();

  /**
   * Creates a protein_index that contains the offset and protein
   * index of the protein in the fasta file.
   */
  ProteinIndex(
    unsigned long int offset, ///< The file location in the database source file
    unsigned int protein_idx ///< The index of the protein in its database.
    );
  
  /**
   * Free a protein index object.
   */
  virtual ~ProteinIndex();

  /**
   * Creates a protein index for the output_file.
   * \returns True if successfully creates a protein index, else
   * false.
   */
  static bool create(
    char* fasta_file ///< input fasta file stream -in
    );

  /**
   * input is the fasta file name which the protein index
   * should have been created.
   * or if creating binary fasta file, is that already on disk?
   *
   *\returns TRUE if protein index or binary fasta file is on disk, else FALSE
   */
  static bool onDisk(
    const char* fasta_file, ///< input fasta file -in
    bool is_binary ///< looking for the binary fasta file or protein index?
    );

};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
