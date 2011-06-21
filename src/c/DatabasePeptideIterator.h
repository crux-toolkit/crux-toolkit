/**
 * \file DatabasePeptideIterator.h 
 * $Revision: 1.27 $
 * \brief Object to iterator over the peptides within a database
 *****************************************************************************/

#ifndef DATABASEPEPTIDEITERATOR_H
#define DATABASEPEPTIDEITERATOR_H

#include "objects.h"
#include "DatabaseProteinIterator.h"
#include "ProteinPeptideIterator.h"
#include "PeptideConstraint.h"

class DatabasePeptideIterator {
 protected:
  DatabaseProteinIterator* database_protein_iterator_; 
    ///< The protein iterator. 
  ProteinPeptideIterator* cur_protein_peptide_iterator_; 
    ///< The peptide iterator for the current protein.
  PeptideConstraint* peptide_constraint_; 
    ///< The constraint for the kind of peptide to iterate over.
  Protein* prior_protein_; 
    ///< the protein that was used before the current working protein
  bool first_passed_; 
    ///< is it ok to convert prior_protein to light?
  PEPTIDE_T* cur_peptide_; ///< the peptide to return by next()
  bool store_all_peptides_; ///< true for search so duplicates are combined
  map<char*, PEPTIDE_T*, cmp_str> peptide_map_; ///< store peptides by sequence
  map<char*, PEPTIDE_T*>::iterator cur_map_position_; ///< next in map to return

  /* Private Functions */
  PEPTIDE_T* nextFromFile();
  bool hasNextFromFile();
  void generateAllPeptides();
  void queueFirstPeptideFromMap();

 public:

  /**
   * Instantiates a new database_peptide_iterator from a database.
   * \returns a DATABASE_PEPTIDE_ITERATOR_T object.
   */
  DatabasePeptideIterator(
    Database* database, ///< the database of interest -in
    PeptideConstraint* peptide_constraint, ///< the peptide_constraint to filter peptides -in
    bool store_all_peptides ///< true: parse all (unique) peptides into map
    );

  /**
   * Frees an allocated database_peptide_iterator object.
   */
  ~DatabasePeptideIterator();

  /**
   * The basic iterator functions.
   * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
   */
  bool hasNext();

  /**
   * \returns The next peptide in the database.
   */
  PEPTIDE_T* next();

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
