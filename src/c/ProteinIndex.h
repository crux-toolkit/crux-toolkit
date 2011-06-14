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
#include "peptide.h"
#include "Protein.h"
#include "carp.h"
#include "PeptideConstraint.h"
#include "sorter.h"


class ProteinIndex {
 protected:
  unsigned long int offset_;  ///< The offset of the protein in the fasta file
  unsigned int protein_idx_;   ///< The protein idx of the protein in the fasta file
 public:


  /**
   * \returns An (empty) protein_index object.
   */
  void init();
  ProteinIndex();

  /**
   * creates a protein_index that contains the offset and protein index of the protein
   * in the fasta file.
   *\returns a new protein_index object
   */
  ProteinIndex(
    unsigned long int offset, ///< The file location in the source file in the database
    unsigned int protein_idx ///< The index of the protein in it's database.
    );
  
  /**
   *
   * free a protein index object
   */
  virtual ~ProteinIndex();

  /**
   * creates a protein index on to the output_file
   * \returns true if successfully creates a protein index, else false
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
    char* fasta_file, ///< input fasta file -in
    bool is_binary ///< are we looking for the binary fasta file? or protein index
    );

};


/**
 * creates a binary fasta file on to the output_file
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
bool create_binary_fasta(
  char* fasta_file  ///< input fasta file -in
  );

/**
 * creates a binary fasta file on to the output_file in currenty directory
 * sets the output file name to the pointer passed in as argument
 * \returns TRUE if successfully creates a binary fasta file, else false
 */
bool create_binary_fasta_in_cur(
  char* fasta_file_w_path, ///< input fasta file with full path -in
  char* fasta_filename, ///< input fasta a file, only filename -in
  char** output_file_name ///< get output filename -out
  );

/**
 * wrapper for create_binary_fasta_file so that two filenames are
 * passed instead of a filename and a filestream.  Eventually should
 * merge to one method
 */
bool create_binary_fasta_here(
  const char* fasta_filename,
  const char* binary_filename
);

/**
 * Heap allocated char*, user must free
 *\returns the binary fasta name which was created from the given fasta file
 */
char* get_binary_fasta_name(
  const char* fasta_file  ///< input fasta file -in                            
  );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
