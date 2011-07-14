/**
 * \file IndexPeptideIterator.h 
 * $Revision: 1.25 $
 * \brief Object for iterating peptides from an index.
 *****************************************************************************/
#ifndef INDEXPEPTIDEITERATOR_H 
#define INDEXPEPTIDEITERATOR_H

#include "Index.h"
#include "IndexFile.h"

class IndexPeptideIterator {
 protected:
  Index* index_; ///< The index object which we are iterating over
  IndexFile* index_files_[MAX_INDEX_FILES]; 
  ///< the index file array that contain information of each index file 
  int total_index_files_; ///< the total count of index_files
  int current_index_file_; ///< the index file open or one to open next 
  FILE* index_file_; ///< The current file stream that we are reading from
  bool has_next_; ///< Is there another peptide?
  PEPTIDE_T* peptide_; ///< the next peptide to return

  /**
   * \brief Parses the "crux_index_map" file that contains the mapping
   * between each crux_index_* file and a mass range. Adds all
   * crux_index_* files that are within the peptide constraint mass
   * range. 
   * \returns true if successfully parses crux_index_map
   */
  bool parseCruxIndexMap();

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
  bool queueNextPeptide();

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

  /**
   * \brief Adds a new index_file object to the index_file.  Checks that
   * the total number of files does not exceed the limit.  Increases the
   * total_index_files count.
   * \returns true if successfully added the new index_file
   */
  bool addNewIndexFile(
    ///< the index_peptide_iterator to add file -out
    char* filename_parsed,  ///< the filename to add -in
    FLOAT_T start_mass,  ///< the start mass of the index file  -in
    FLOAT_T range  ///< the mass range of the index file  -in
    );

 public:

  /**
   * Instantiates a new peptide_iterator from an index, which returns peptides
   * that obey peptide constraint. At first will only accept constraints
   * that will require reading in one file (e.g a 1m/z range). 
   */
  IndexPeptideIterator(
    Index* index ///< The index -in
    );

  /**
   * Frees an allocated index_peptide_iterator object.
   */
  ~IndexPeptideIterator();

  /**
   * The basic iterator functions.
   * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
   */
  bool hasNext();

  /**
   * \returns The next peptide in the index.
   */
  PEPTIDE_T* next();

  Index* getIndex();



};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
