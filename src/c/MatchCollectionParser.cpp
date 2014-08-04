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
#include "SQTReader.h"
#include "MzIdentMLReader.h"
#include "Protein.h"
#include "PostProcessProtein.h"
#include <sys/types.h>
#include <sys/stat.h>

using namespace Crux;

MatchCollectionParser::MatchCollectionParser() {
    database_ = NULL;
    decoy_database_= NULL;
}

MatchCollectionParser::~MatchCollectionParser() {
    Database::freeDatabase(database_);
    Database::freeDatabase(decoy_database_);
}


/**
   * Creates database object(s) from fasta or index file
   */
void MatchCollectionParser::loadDatabase(
  const char* fasta_file, ///< fasta or index path -in  
  Database*& database, ///< resulting database -out
  Database*& decoy_database ///< resulting decoy database -out
  ) {

  if (fasta_file == NULL || string(fasta_file) == string("__NULL_STR")) {
    carp(CARP_DEBUG, "no database provided");
    database = new Database();
    database -> setIsParsed(true);
    decoy_database = new Database();
    database -> setIsParsed(true);
  } else {
    bool use_index = is_directory(fasta_file);
    // get binary fasta file name with path to crux directory 
    if (use_index == true){ 
      //We aren't supporting search-for-matches anymore, so we won't be
      //building any indices.  We could use tide-indices for this, but
      //comet doesn't use an index.  In the future we can support tide
      //indices in the post-processing tools if users want that.
      carp(CARP_FATAL, "Parsing binary index from tide-index is not"
        " supported.\nPlease use the original fasta file instead.");
    } else {
      database = new Database(fasta_file, false);// not memmapped
      database->transformTextToMemmap(".", true);// is temp
      decoy_database = new Database();
    }
    database->parse();
  }
}


/**
 * \returns a Protein object.  This function first tries to
 * find the protein by id in the target or decoy database.
 * if unsucessful, we then determine whether the protein
 * is a decoy or not by looking for the "decoy-prefix" in
 * id.  When then create a PostProcessProtein and add it
 * to the corresponding target or decoy database.
 */
Protein* MatchCollectionParser::getProtein(
  Database* database, ///< target database -in
  Database* decoy_database, ///< decoy database -in
  string& protein_id, ///< id of protein to find -in
  bool& is_decoy ///< is protein a decoy? -out
  ) {

  //is it in the target database?
  Protein* protein = database->getProteinByIdString(protein_id.c_str());
  if (protein != NULL) { 
    is_decoy = false;
    return protein; 
  }

  //is it in the decoy database?
  if (decoy_database != NULL) {
    protein = decoy_database->getProteinByIdString(protein_id.c_str());
    if (protein != NULL) {
      is_decoy = true;
      return protein;
    }
  }

  //try creating it and adding it to the database as a postprocess protein
  carp(CARP_DEBUG, "Creating new protein for %s",protein_id.c_str());
  protein = new PostProcessProtein();
  protein->setId(protein_id.c_str());
  string decoy_prefix = get_string_parameter_pointer("decoy-prefix");
  if (protein_id.find(decoy_prefix) != string::npos) {
    carp(CARP_DEBUG, "adding to decoy database");
    is_decoy = true;
    decoy_database->addProtein(protein);
    
  } else {
    carp(CARP_DEBUG, "adding to target database");
    is_decoy = false;
    database->addProtein(protein);
  }
  return protein;
}

Protein* MatchCollectionParser::getProtein(
  Database* database, ///< target database -in
  Database* decoy_database, ///< decoy database -in
  string& protein_id, ///< id of protein to find -in
  string& sequence, ///< sequence of the protein -in
  bool& is_decoy ///< is protein a decoy? -out
  ) {

  //is it in the target database?
  Protein* protein = database->getProteinByIdString(protein_id.c_str());
  if (protein != NULL) { 
    is_decoy = false;
    return protein; 
  }

  //is it in the decoy database?
  if (decoy_database != NULL) {
    protein = decoy_database->getProteinByIdString(protein_id.c_str());
    if (protein != NULL) {
      is_decoy = true;
      return protein;
    }
  }

  //try creating it and adding it to the database as a postprocess protein
  carp(CARP_DEBUG, "Creating new protein for %s",protein_id.c_str());
  carp(CARP_DEBUG, "Sequence :%s",sequence.c_str());
  protein = new Protein();
  protein->setId(protein_id.c_str());
  protein->setSequence(sequence.c_str());
  protein->setLength(sequence.length());
  

  string decoy_prefix = get_string_parameter_pointer("decoy-prefix");
  if (protein_id.find(decoy_prefix) != string::npos) {
    is_decoy = true;
    decoy_database->addProtein(protein);
  } else {
    is_decoy = false;
    database->addProtein(protein);
  }
  return protein;
}

/**
 * \returns a MatchCollection object using the file and protein database
 */
MatchCollection* MatchCollectionParser::create(
  const char* match_path, ///< path to the file of matches 
  const char* fasta_path ///< path to the protein database
  ) {

  carp(CARP_DEBUG, "match path:%s", match_path);
  if (fasta_path) {
    carp(CARP_DEBUG, "fasta path:%s", fasta_path);
  }
  struct stat stat_buff ; 
  stat(match_path, &stat_buff);
  if(stat(match_path, &stat_buff)!=0){
    carp(CARP_FATAL, "The file %s does not exist. \n", match_path);
  }
  
  if (database_ == NULL || decoy_database_ == NULL) {
    loadDatabase(fasta_path, database_, decoy_database_);
  }
  MatchCollection* collection = NULL;
  
  if (S_ISDIR(stat_buff.st_mode)){
    carp(CARP_FATAL, "Internal error");
  }else if( has_extension(match_path, ".txt")){
    collection = MatchFileReader::parse(match_path, database_, decoy_database_);
  } else if( has_extension(match_path, ".xml")) {
    collection = PepXMLReader::parse(match_path, database_, decoy_database_);
  } else if (has_extension(match_path, ".sqt")) {
    collection = SQTReader::parse(match_path, database_, decoy_database_);
  } else if( has_extension(match_path, ".mzid")) {
    collection = MzIdentMLReader::parse(match_path, database_, decoy_database_);
  } else {
    collection = MatchFileReader::parse(match_path, database_, decoy_database_);
  }
  collection->setFilePath(match_path);

  return collection;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
