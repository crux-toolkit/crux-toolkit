/**
 * \file MatchCollectionIterator.h 
 * $Revision: 1.38 $
 * \brief An object that iterates over the MatchCollection objects in
 * the specified directory of serialized MatchCollections 
 */
#ifndef MATCHCOLLECTIONITERATOR_H
#define MATCHCOLLECTIONITERATOR_H

#include "MatchCollection.h"

/******************************
 * match_collection_iterator
 ******************************/
class MatchCollectionIterator {
 protected:
  DIR* working_directory_; 
  ///< the working directory for the iterator to find match_collections
  char* directory_name_; ///< the directory name in char
  Database* database_; ///< database of proteins in the results files
  Database* decoy_database_; ///< optional database of decoy proteins in files
  int number_collections_; 
  ///< the total number of match_collections in the directory (target+decoy)
  int collection_idx_;  ///< the index of the current collection to return
  MatchCollection* match_collection_; ///< the match collection to return
  bool has_another_collection_; ///< has another match_collection to return?
  vector<bool>* cols_in_file_; ///< which columns were in the target file

  /**
   * \brief Finds the next match_collection in directory and prepares
   * the iterator to hand it off when 'next' called.
   *
   * When no more match_collections (i.e. psm files) are available, set
   * match_collection_iterator->has_another_collection to false
   * \returns void
   */
  void setup();

 public:

  /**
   * Create a match_collection iterator from a directory of serialized files
   * Only hadles up to one target and three decoy sets per folder
   *\returns match_collection iterator instantiated from a result folder
   */
  MatchCollectionIterator(
    const char* output_file_directory, ///< the directory path where the PSM output files are located -in
    const char* fasta_file, ///< The name of the file (in fasta format) from which to retrieve proteins and peptides for match_collections. -in
    int* decoy_count
    );

  /**
   * free match_collection_iterator
   */
  virtual ~MatchCollectionIterator();


  /**
   *\returns true, if there's another match_collection to return, else return false
   */
  bool hasNext();
  
  /**
   * returns the next match collection object and sets up fro the next iteration
   *\returns the next match collection object
   */
  MatchCollection* next();

  /**
   *\returns the database
   */
  Database* getDatabase();
    
  /**
   *\returns the decoy database
   */
  Database* getDecoyDatabase();

  /**
   *\returns the total number of match_collections to return
   */
  int getNumberCollections();
  
  /**
   * \brief Get the name of the directory the match_collection_iterator
   * is working in.
   * \returns A const pointer to the directory name.
   */
  const char* getDirectoryName();
  
  /**
   * \brief Get the working directory
   */
  DIR* getWorkingDirectory();



  vector<bool>& getColsInFile();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
