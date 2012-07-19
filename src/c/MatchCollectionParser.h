/**
 * \file MatchCollectionParser.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 14 June 2011
 * \brief Return a MatchCollection object from the appropriate
 * file format.
 */
#ifndef MATCHCOLLECTIONPARSER_H
#define MATCHCOLLECTIONPARSER_H

#include "MatchCollection.h"

/**
 * Instantiates a MatchCollection based on the extension of the
 * given file.
 */
class MatchCollectionParser {

 protected:

 public:

  /**
   * \returns a MatchCollection object using the file and protein database
   */
  static MatchCollection* create(
    const char* match_path, ///< path to the file of matches
    const char* fasta_path  ///< path to the protein database
  );


  /**
   * Creates database object(s) from fasta or index file
   */
  static void loadDatabase(
    const char* fasta_file, ///< fasta or index path -in
    Database*& database,  ///< resulting database -out
    Database*& decoy_database ///< resulting decoy database -out
  );

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif 
