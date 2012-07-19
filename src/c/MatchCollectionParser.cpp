/**
 * \file MatchCollectionParser.cpp 
 * AUTHOR Sean McIlwain
 * CREATE DATE: 14 June 2011
 * \brief Return a MatchCollection object of the appropriate
 * derived class.
 */
#include "parameter.h"
#include "MatchCollectionParser.h"
#include "MatchFileReader.h"
#include "PepXMLReader.h"
#include <sys/types.h>
#include <sys/stat.h>

/**
   * Creates database object(s) from fasta or index file
   */
void MatchCollectionParser::loadDatabase(
  const char* fasta_file, ///< fasta or index path -in  
  Database*& database, ///< resulting database -out
  Database*& decoy_database ///< resulting decoy database -out
  ) {

  bool use_index = is_directory(fasta_file);
  // get binary fasta file name with path to crux directory 
  if (use_index == true){ 
    char* binary_fasta = Index::getBinaryFastaName(fasta_file);
    database = new Database(binary_fasta, true);// is memmapped
    free(binary_fasta);
    binary_fasta = Index::getDecoyBinaryFastaName(fasta_file);
    if( binary_fasta != NULL ){
      decoy_database = new Database(binary_fasta, true);// is memmapped
      decoy_database->parse();
      free(binary_fasta);
    }
  } else {
    database = new Database(fasta_file, false);// not memmapped
    database->transformTextToMemmap(".", true);// is temp
    decoy_database = NULL;
  }
  database->parse();

}


/**
 * \returns a MatchCollection object using the file and protein database
 */
MatchCollection* MatchCollectionParser::create(
  const char* match_path, ///< path to the file of matches 
  const char* fasta_path ///< path to the protein database
  ) {

  struct stat stat_buff ; 
  stat(match_path, &stat_buff);
  if(stat(match_path, &stat_buff)!=0){
    carp(CARP_FATAL, "The file %s does not exist. \n", match_path);
  }

  Database* database;
  Database* decoy_database;
  
  loadDatabase(fasta_path, database, decoy_database);

  MatchCollection* collection = NULL;

  if (S_ISDIR(stat_buff.st_mode)){
    carp(CARP_FATAL, "Internal error");
  }else if( has_extension(match_path, ".txt")){
    collection = MatchFileReader::parse(match_path, database, decoy_database);
  } else if( has_extension(match_path, ".xml")) {
    collection = PepXMLReader::parse(match_path, database, decoy_database);
  } else {
    collection = MatchFileReader::parse(match_path, database, decoy_database);
  }
   
  return collection;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
