/**
 * \file IndexPeptideIterator.h 
 * $Revision: 1.25 $
 * \brief Object for iterating peptides from an index.
 *****************************************************************************/
#ifndef INDEXPEPTIDEITERATOR_H 
#define INDEXPEPTIDEITERATOR_H

#include "Index.h"
#include "IndexFile.h"
#include "PeptideIterator.h"

#include <vector>

class IndexPeptideIterator : public PeptideIterator {
 protected:
  Index* index_; ///< The index object which we are iterating over
  std::vector<IndexFile*> index_files_;  
  ///< the index file vector that contains information of each index file 
  size_t current_index_file_; ///< the index file open or one to open next 
  FILE* index_file_; ///< The current file stream that we are reading from
  bool is_decoy_;    ///< return decoy or target peptides

  PeptideConstraint* constraint_;

  /**
   * \brief Find the next peptide for the iterator to return.
   *
   * Called in the process of initializing a new iterator and by
   * next().  Looks in the current index_file for peptides matching the
   * constraint.  If none found, checks remaining index_files.  If no
   * constraint-satisfying peptides are found in any of the remaining
   * index files, next_peptide is set to NULL and has_next is set to false.
   * \returns true if no errors were encountered while reading files
   * (even if there is no peptide to return).
   */
  virtual bool queueNextPeptide();

  /**
   * \brief Prepare an index peptide iterator to have its index files
   * searched.  Changes the state of the iterator if its index_file
   * is NULL and the current file is not the last in the list.  The
   * routine that searches for peptides in a file closes any files it
   * has read to the end and sets the pointer to NULL.
   * \returns true if there is a file ready to be read or false if no
   * more files remain.
   */
  bool findNextIndexFile();

  /**
   * \brief Search for a peptide matching the peptide constraint in the
   * current index file of the iterator.  If the iterator has no open
   * index file, returns false.  If the current file reaches the end
   * without finding a suitable peptide, the file is closed, file handle
   * set to NULL, and the number of the current file is increemnted.
   *
   * \returns true if the iterator is ready to return a new peptide,
   * false if no peptide meeting the constraint could abe found in the
   * current file. 
   */
  bool findPeptideInCurrentIndexFile();

  /**
   * \brief Find the next peptide in the file that meets the peptide
   * constraint, read it from file, and add set up the iterator to
   * return it.
   * \returns true if successfully finds and parses a peptide that meets
   * the constraint.
   */
  bool fastForwardIndexFile();

 public:

  /**
   * Instantiates a new peptide_iterator from an index, which returns peptides
   * that obey peptide constraint. 
   */
  IndexPeptideIterator(
    Index* index, ///< The index -in
    PeptideConstraint* constraint, ///< the constraint -in
    bool is_decoy = false ///< return target or decoy peptides
    );

  /**
   * Frees an allocated index_peptide_iterator object.
   */
  ~IndexPeptideIterator();

  Index* getIndex();

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
