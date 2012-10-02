/**
 * \file IndexFilteredPeptideIterator.h 
 * $Revision: 1.25 $
 * \brief Object for iterating peptides from an index.  Unlike a
 * standard index iterator, it filters by digestion type
 *****************************************************************************/
#ifndef INDEX_FILTERED_PEPTIDE_ITERATOR_H 
#define INDEX_FILTERED_PEPTIDE_ITERATOR_H 

#include "IndexPeptideIterator.h"

/**
 * \class IndexFilteredPeptideIterator
 * \brief Returns peptides from an index that have been further
 * filtered by digestion type.
 */
class IndexFilteredPeptideIterator : public IndexPeptideIterator {
 protected:
  bool has_next_;      ///< Are there more peptides to return?
  Crux::Peptide* peptide_; ///< The next peptide to return

  /**
   * Sets up the iterator.
   * \returns true if successfully sets up the iterator, else false.
   */
  bool setup();

 public:
  /**
   * Instantiates a new index_filtered_peptide_iterator from a index.
   * \returns a new heap allocated index_filtered_peptide_iterator object
   */
  IndexFilteredPeptideIterator(
    Index* index ///< The index object which we are iterating over -in
  );

  /**
   * Frees an allocated index_filtered_peptide_iterator object.
   */
  ~IndexFilteredPeptideIterator();

  /**
   *  The basic iterator functions.
   * \returns The next peptide in the index.
   */
  Crux::Peptide* next();

  /**
   * The basic iterator functions.
   * Does the iterator have more peptides to return.
   *\returns true if there are additional peptides, false if not.
   */
  bool hasNext();


};

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Crux::Peptide* void_index_filtered_peptide_iterator_next(
  void* index_filtered_peptide_iterator ///< the iterator to initialize -in
  );

/**
 * The basic iterator functions.
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
bool void_index_filtered_peptide_iterator_has_next(
  void* index_filtered_peptide_iterator ///< the iterator to initialize -in
  );

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void void_free_index_filtered_peptide_iterator(
  void* index_filtered_peptide_iterator ///< the iterator to initialize -in
  );


#endif //INDEX_FILTERED_PEPTIDE_ITERATOR_H 
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
