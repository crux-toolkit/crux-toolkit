/**
 * \file IndexMap.h
 * $Revision: 1.00 $
 * \brief Object representing the index map built for an index
 */

#ifndef INDEXMAP_H_
#define INDEXMAP_H_

#include <vector>
#include "objects.h"

#include "IndexFile.h"

class IndexMap {
 protected:
  Index* index_;
  bool is_decoy_;

  std::vector<IndexFile*> index_files_; ///< Should have all the index files in a index

  /**
   * \brief Parses the "crux_index_map" file that contains the mapping
   * between each crux_index_* file and a mass range. Adds all
   * crux_index_* files that are within the peptide constraint mass
   * range. 
   * \returns true if successfully parses crux_index_map
   */
  bool parse(); 

  /**
   * \brief Adds a new index_file object to the index_file.  Checks that
   * the total number of files does not exceed the limit.  Increases the
   * total_index_files count.
   * \returns true if successfully added the new index_file
   */
  bool addNewIndexFile(
    const char* filename_parsed,  ///< the filename to add -in
    FLOAT_T start_mass,  ///< the start mass of the index file  -in
    FLOAT_T range  ///< the mass range of the index file  -in
    );

 public:
  
  /**
   * Loads the index map into memory
   */
  IndexMap(
    Index* index, ///< Index where the map should be generated from
    bool is_decoy ///< Are we using the decoy or target index?
  );

  /**
   * Default destructor
   */
  virtual ~IndexMap();

  
  /**
   * populates the constraint_index_files with the index files that match
   * the passed in constraint
   */
  void getIndexFiles(
    PeptideConstraint* constraint, ///< the constraint for peptides to select -in 
    std::vector<IndexFile*>& constraint_index_files ///< a list of all of the indices passing the constraint -out
  );

};


#endif
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
