/*************************************************************************//**
 * \file Database.cpp
 * \brief Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#include <fcntl.h>
#ifndef _MSC_VER
#include <sys/mman.h>
#include <unistd.h>
#endif
#include "utils.h"
#include "crux-utils.h"
#include "Peptide.h"
#include "Protein.h"
#include "ProteinPeptideIterator.h"
#include "Database.h"
#include "hash.h"
#include "carp.h"
#include "objects.h"
#include "PeptideConstraint.h"
#include "ProteinIndex.h"
#include "ProteinIndexIterator.h"

#include "DatabaseProteinIterator.h"
#include "DatabasePeptideIterator.h"

#include <map>
#include <vector>
#include <iostream>

using namespace std;
using namespace Crux; 
const string Database::binary_suffix = "-binary-fasta";
const string Database::decoy_binary_suffix = "-binary-fasta-decoy";
const string Database::decoy_fasta_suffix = "-random.fasta";

/**
 * intializes a database object
 */
void Database::init(){
  file_ = NULL;
  is_parsed_ = false;
  size_ = 0; 
  use_light_protein_ = false; 
  is_memmap_ = false;
  data_address_ = NULL;
  pointer_count_ = 1;
  file_size_ = 0;
  is_hashed_ = false;
  proteins_ = new vector<Protein*>();
  protein_map_ = new map<char*, Protein*, cmp_str>();
  decoys_ = NO_DECOYS;
  binary_is_temp_ = false;
}

/**
 * \returns An (empty) database object.
 */
Database::Database() {
  init();
}


/**
 * \returns A new database object.
 */
Database::Database(
  const char*         filename, ///< The file from which to parse the database. 
  ///< either text fasta file or binary fasta file -in
  bool is_memmap, ///< are we using a memory mapped binary fasta file? 
  ///< If so, all proteins are memory mapped -in
  DECOY_TYPE_T decoys ///< is this a decoy database
  )         
{
  carp(CARP_DEBUG, "Creating new database from '%s'", filename);
  init();
  is_memmap_ = is_memmap;
  if( is_memmap_ ){
    binary_filename_ = filename;
  } else {
    fasta_filename_ = filename;
  }
  decoys_ = decoys;
}  

/**
 * Frees an allocated protein object.
 */
void Database::freeDatabase(
  Database* database ///< An allocated database -in
  )
{
  if(database == NULL){
    return;
  }

  // decrement database pointer counter
  --database->pointer_count_;
  carp(CARP_DETAILED_DEBUG, "Database pointer count %i",  
      database->pointer_count_);

  // DEBUG show the databse pointer count
  // printf("Free: After free: %s: %d\n", database->pointer_count);

  // only free up memory when remaining pointers are from the proteins
  if((size_t)database->pointer_count_ > database->proteins_->size() ){ 
    return;
  }

  delete database;
}

Database::~Database() {

  // only free proteins if been parsed and file has been opened
  if(is_parsed_){
    carp(CARP_DEBUG, "Freeing database.");
    
    // free each protein in the array
    unsigned int protein_idx;
    for(protein_idx=0; protein_idx < proteins_->size(); ++protein_idx){
      delete ((*proteins_)[protein_idx]);
    }
    delete proteins_;
    delete protein_map_; // contents already deleted
    
    // free memory mapped binary file from memory
    if(is_memmap_){
      // un map the memory!!
#ifdef _MSC_VER
      stub_unmmap(&unmap_info_); 
#else
      if(munmap(data_address_, file_size_) != 0){
        carp(CARP_ERROR, "failed to unmap the memory of binary fasta file");
      }
#endif
    }
    // not memory mapped
    else if (file_ != NULL) {
      // close file handle
      carp(CARP_DEBUG, "Closing database filehandle");
      fclose(file_);
    }
  }

  if( binary_is_temp_ && !binary_filename_.empty() ){
    carp(CARP_DEBUG, "Deleting temp binary fasta %s.", 
         binary_filename_.c_str());
    remove(binary_filename_.c_str());
  }
}

/**
 * Prints a database object to file.
 */
void Database::print(
  FILE* file    ///< output file stream -out             
  )
{
  Protein* protein = NULL;

  if( is_memmap_ ){
    fprintf(file, "filename:%s\n", binary_filename_.c_str());
  } else {
    fprintf(file, "filename:%s\n", fasta_filename_.c_str());
  }
  fprintf(file, "is_parsed:");
  
  // has the database been parsed?
  if(is_parsed_){
    fprintf(file, "true\n");
    DatabaseProteinIterator* iterator
      = new DatabaseProteinIterator(this);
 
    while(iterator->hasNext()){
      protein = iterator->next();
      // if the database uses light/heavy functionality
      if(use_light_protein_){
        protein->toHeavy();
      }
      protein->print(stdout);
  
      // if the database uses light/heavy functionality
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      if(database->use_light_protein){
        protein_to_light(protein);
      }
      */
    }
    delete iterator;
  }
  else{
    fprintf(file, "false\n");
  }
}

void Database::addProtein(
  Protein* protein
  ) {

  protein->setDatabase(this);
      
  // add protein to database
  proteins_->push_back(protein);
  
  protein->setProteinIdx(proteins_->size()-1);

  if (is_hashed_) {
    char* id = protein->getIdPointer();
    protein_map_->insert(make_pair(id, protein));
  }
}

/**
 * Parses a database from the text based fasta file in the filename
 * member variable
 * reads in all proteins in the fasta file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 * IF using light_protein functionality will not read in the sequence or id.
 * \returns true if success. false if failure.
 */
bool Database::parseTextFasta()
{
  unsigned long working_index;
  FILE* file = NULL;
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  Protein* new_protein;
  unsigned int protein_idx;

  carp(CARP_DEBUG, "Parsing text fasta file '%s'", fasta_filename_.c_str());
  // check if already parsed
  if(is_parsed_){
    return true;
  }
  
  // open file and 
  file = fopen(fasta_filename_.c_str(), "rb");
  
  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_ERROR, "Failed to open fasta file %s", fasta_filename_.c_str());
    return false;
  }
  
  // check if use light protein and parse those proteins from protein index
  if(use_light_protein_ && ProteinIndex::onDisk(fasta_filename_.c_str(), false)){
    // let the user know that protein index file is being used
    carp(CARP_INFO, "using protein index file");

    // create a protein index iterator
    ProteinIndexIterator* protein_index_iterator =
      new ProteinIndexIterator(fasta_filename_.c_str());

    // iterate over all proteins in protein index
    while(protein_index_iterator->hasNext()){
      // check if there's space for more proteins
      new_protein = protein_index_iterator->next();
      new_protein->setDatabase(this);
      
      // add protein to database
      proteins_->push_back(new_protein);
    }
    // job well done..free iterator
    delete protein_index_iterator;
  }
  else{  
    working_index = ftell(file);
    // check each line until reach '>' line
    while((line_length =  getline(&new_line, &buf_length, file)) != -1){
      if(new_line[0] == '>'){
        // the new protein to be added
        new_protein = new Protein();
        
        // do not parse the protein sequence if using light/heavy functionality
        if(use_light_protein_){
          // set light and offset
          new_protein->setOffset(working_index);
          new_protein->setIsLight(true);
        }
        else{
          // rewind to the beginning of the protein to include ">" line
          fseek(file, working_index, SEEK_SET);
          
          // failed to parse the protein from fasta file
          // protein offset is set in the parse_protein_fasta_file method
          if(!new_protein->parseProteinFastaFile(file)){
            fclose(file);
            delete new_protein;
            for(protein_idx=0;protein_idx<proteins_->size();protein_idx++){
              delete (proteins_->at(protein_idx));
            }
            proteins_->clear();
            carp(CARP_ERROR, "Failed to parse fasta file");
            return false;
          }
          new_protein->setIsLight(false);
        }
        // add protein to database
        proteins_->push_back(new_protein);
        // set protein index, database
        new_protein->setProteinIdx(proteins_->size()-1);
        new_protein->setDatabase(this);
      }
      working_index = ftell(file);
    }
    free(new_line);
  }
  
  // yes the database is paresed now..!!
  is_parsed_ = true;
  file_ = file;
  return true;
}

/**
 * memory maps the binary fasta file for the database
 *\return true if successfully memory map binary fasta file, else false
 */
bool Database::memoryMap(
  int file_d  ///<  file descriptor -in
  )
{
  struct stat file_info;
  
  // get information of the binary fasta file
  if (stat(binary_filename_.c_str(), &file_info) == -1) {
    carp(CARP_ERROR,
         "Failed to retrieve information of binary fasta file: %s",
         binary_filename_.c_str());
    return false;
  }
  
  // set size of the binary fasta file in database
  // this is used later to know how much to unmap
  file_size_ = file_info.st_size;
  
  // memory map the entire binary fasta file!
#ifdef _MSC_VER
  data_address_ = stub_mmap(binary_filename_.c_str(), &unmap_info_);

#else
  data_address_ = mmap((caddr_t)0, file_info.st_size, 
                       PROT_READ, MAP_PRIVATE /*MAP_SHARED*/, file_d, 0);

  // check if memory mapping has succeeded
  if ((caddr_t)(data_address_) == (caddr_t)(-1)){
    carp(CARP_ERROR, "Failed to use mmap function for binary fasta file: %s", 
         binary_filename_.c_str());
    return false;
  }
#endif
  
  return true;
}

/**
 * Assumes that there is a 1 at the very end after all the proteins in
 * binary file.
 * \returns True successfully populates the proteins from memory mapped
 * binary fasta file, else false .
 */
bool Database::populateProteinsFromMemmap()
{
  Protein* new_protein;
  unsigned int protein_idx = 0;
  char* data = (char*)data_address_;
  
  // parse proteins until the end of list
  while((int)data[0] != 1){
    
    // the new protein to be added
    new_protein = new Protein();
    
    // parse protein from memory map
    if(!new_protein->parseProteinBinaryMemmap(&data)){
      // failed to parse the protein from memmap
      // free all proteins, and return false
      delete new_protein;
      for(; protein_idx < proteins_->size(); ++protein_idx){
        delete proteins_->at(protein_idx);
      }
      proteins_->clear();
      carp(CARP_ERROR, "Failed to parse fasta file");
      return false;
    }
    new_protein->setIsLight(false);
    
    // add protein to database
    proteins_->push_back(new_protein);
    // set protein index, database
    new_protein->setProteinIdx(proteins_->size()-1);
    new_protein->setDatabase(this);
  }

  return true;
}

/**
 * \brief Parses a database from the binary fasta file in the filename
 * member variable.
 *
 * Memory maps the binary fasta file into memory. The protein
 * sequences are not copied, but just pointed to the memory mapped
 * location. 
 * \returns true if success. false if failure.
 */
bool Database::parseMemmapBinary()
{
  int file_d = -1;
  carp(CARP_DEBUG, "Parsing binary fasta file '%s'", binary_filename_.c_str());
 
  // check if already parsed
  if(is_parsed_){
    return true;
  }
  
  // open file and 
  file_d = open(binary_filename_.c_str(), O_RDONLY);
  
  // check if succesfully opened file
  if(file_d == -1){
    carp(CARP_FATAL, "Failed to open file to parse database");
    return false;
  }

  // FIXME, if what to use some light protein for binary file change here...
  // check if user request light protein
  // When using a binary file in memory map, cannot use light protein
  // change to false on light protein useage
  if(use_light_protein_){
    carp(CARP_WARNING, 
         "memory mapping does not support light protein,changing settings to use heavy protein");
    use_light_protein_ = false;;
  }

  // memory map the binary fasta file into memory
  if(!memoryMap(file_d)){
    carp(CARP_ERROR, "Failed to memory map binary fasta file into memory");
    return false;
  }

  // populate the proteins from the memory mapped fasta file
  if(!populateProteinsFromMemmap()){
    carp(CARP_ERROR, "Failed to populate the proteins from memory mapped fasta file");
    return false;
  }
   
  // yes the database is paresed now..!!
  is_parsed_ = true;
  return true;
}


/**
 * Parses a database from the file in the filename member variable
 * The is_memmap field in the database struct determines whether the
 * input file is a binary fasta file or normal text fasta file.
 *
 * IF is_memmap is true, memory maps the entire binary fasta file into memory
 * and then creates protein objects that point to the memory mapped binary file
 *
 * IF is_memmap is false, uses the traditional text fasta file which
 * it parses out the various peptides for each protein. Only when
 * using text fasta file can you use light/heavy protein, in which 
 * if using light_protein functionality will not read in the sequence
 * or id. Will parse sequence if protein  
 * is needed, lazy parsing.
 *
 * For Both cases, reads in all proteins in file and creates a protein object
 * and adds them to the database protein array
 * total proteins in fasta file must not exceed MAX_PROTEIN constant
 *
 * \returns true if success. false if failure.
 */
bool Database::parse()
{
  // should we parse the database using memory mapped binary fasta file?
  if(is_memmap_){
    return parseMemmapBinary();   
  }
  else{ // parse database from normal text fasta file, no memory mapping!
    return parseTextFasta();
  }
  
  // succeeded in parsing database
  return true;
}


/**
 * \brief Changes a database from one that reads from a fasta file to
 * one that reads from a binary/memmory mapped protein file.
 *
 * If database already has binary source (i.e. is_memmap == true), 
 * returns true.  
 * Opens the fasta file pointed to by filename for reading.  Creates an
 * output file with the name given.  Reads in each protein from the
 * text file and serializes it to the output file.  Closes both files.
 * Changes filename to point to new output file and sets is_memmap to
 * true. Parses the database.
 * \returns true if all processes succeed, else false.
 */
bool Database::transformTextToMemmap(
  const char* output_dir,
  bool binary_is_temp
  ){

  createBinaryFasta(output_dir, binary_is_temp);

  // set is_memmap to true
  setMemmap(true);

  // parse the binary fasta
  bool success = parse();

  return success;
}

/**
 * Using the fasta file the Database was instantiated with, write a
 * binary protein file in the given directory to use for memory
 * mapping.  If is_temp, delete the file on destruction.  Warns if
 * Database was not opened with a text file.  If the database is to
 * contain decoy proteins, randomizes each protein before
 * serializing.  Also prints a new fasta file of the decoy proteins in
 * the same directory as the binary file.
 */
void Database::createBinaryFasta(const char* directory, bool is_temp){
  binary_is_temp_ = is_temp;

  if( fasta_filename_.empty() ){
    carp(CARP_WARNING, "No fasta file to transform to binary.");
    return;
  }

  // get the name of the new binary file
  const char* binary_suffix = Database::binary_suffix.c_str();
  if( decoys_ != NO_DECOYS ){
    binary_suffix = Database::decoy_binary_suffix.c_str();
  }
  char* tmp_name = generate_name_path( fasta_filename_.c_str(), ".fasta",
                                       binary_suffix, directory);
  binary_filename_ = tmp_name;
  free(tmp_name);
  carp(CARP_DEBUG, "Transforming text file '%s' to binary file '%s'",
       fasta_filename_.c_str(), binary_filename_.c_str());

  // open output file
  FILE* output_file = fopen(binary_filename_.c_str(), "wb");
  if( output_file == NULL ){
    carp(CARP_FATAL, "Could not open binary protein file %s", 
         binary_filename_.c_str());
  }
  // also open a fasta file if this is a decoy database
  FILE* output_fasta = NULL;
  if( decoys_ != NO_DECOYS){
    vector<const char*> suffixes;
    suffixes.push_back(".fasta");
    suffixes.push_back(".fa");
    suffixes.push_back(".fsa");
    char* fasta_output_name = generate_name_path( fasta_filename_.c_str(),
                                                  suffixes, 
                                                  decoy_fasta_suffix.c_str(),
                                                  directory);
    output_fasta = fopen(fasta_output_name, "wb");
    if( output_fasta == NULL ){
      carp(CARP_FATAL, "Could not open new fasta file %s for decoy proteins.",
           output_fasta);
    }
    free(fasta_output_name);
  }


  // open input file
  FILE* input_file = fopen(fasta_filename_.c_str(), "rb");
  if( input_file == NULL ){
    carp(CARP_FATAL, "Could not open fasta file %s", fasta_filename_.c_str());
  }

  // for file reading
  char* new_line = NULL;
  int line_length = 0;
  size_t buf_length = 0;
  unsigned int mod_me = 1000;
  unsigned int protein_idx = 0;

  // read through the fasta and at each line beginning with >, parse a protein
  unsigned long working_index = ftell(input_file);
  while((line_length =  getline(&new_line, &buf_length, input_file)) != -1){
    if(new_line[0] == '>'){
      // the new protein to be serialize
      Protein* new_protein = new Protein();
      
      // rewind to the begining of the protein to include ">" line
      fseek(input_file, working_index, SEEK_SET);
      // protein offset is set in the parse_protein_fasta_file method
      if(!new_protein->parseProteinFastaFile(input_file)){
        fclose(input_file);
        delete new_protein;
        carp(CARP_FATAL, "Failed to parse fasta file");
      }
      new_protein->setIsLight(false);

      if( decoys_ != NO_DECOYS ){
        new_protein->shuffle(decoys_);
        new_protein->print(output_fasta);
      }
      
      // serialize protein as binary to output file
      new_protein->serialize(output_file);

      // update protein count
      ++protein_idx;

      // free this protein
      delete new_protein;
    } 

    // print status
    if(protein_idx % mod_me == 0){
      if((protein_idx / 10) == mod_me){
        mod_me *= 10;
      }
      carp(CARP_INFO, "Reached protein %d", protein_idx);
    }

    working_index = ftell(input_file);
  } // next line

  // write the end character to binary fasta file
  char term_char = 1;  // use 1 and not '*' as the terminal
                       // character for the file b/c id length is 
                       // stored in same field and id len == 42
                       // is the smae as '*'
  fwrite(&term_char, sizeof(char), 1, output_file);

  // print final status
  carp(CARP_INFO, "Total proteins found: %d", protein_idx);
  
  free(new_line);
  fclose(input_file);
  fclose(output_file);
}

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 *sets true,false whether the database uses memory mapped
 */
void Database::setMemmap(
  bool is_memmap  ///< is the database memory mapped?
  )
{
  is_memmap_ = is_memmap;
}

/**
 *\returns the filename of the database
 * returns a heap allocated new copy of the filename
 * user must free the return filename
 */
char* Database::getFilename()
{
  if( is_memmap_ ){
    return my_copy_string(binary_filename_.c_str());
  }
  return my_copy_string(fasta_filename_.c_str());
}


/**
 *\returns the pointer to the filename of the database
 * user must not free or change the filename
 */
const char* Database::getFilenamePointer()
{
  if( is_memmap_ ){
    return binary_filename_.c_str();
  }
  return fasta_filename_.c_str();
}

/**
 *\returns true|false whether the database has been parsed?
 */
bool Database::getIsParsed()
{
  return is_parsed_;
}

void Database::setIsParsed(
  bool is_parsed
  ) {
  is_parsed_ = is_parsed;
}

/**
 * \returns The type of shuffling used on the proteins in this database
 */
DECOY_TYPE_T Database::getDecoyType(){
  return decoys_;
}


/**
 * sets the use_light_protein of the database
 */
void Database::setUseLightProtein(
  bool use ///< should I use the light/heavy functionality?
  )
{
  use_light_protein_ = use;
}

/**
 *\returns true|false whether the database uses light/heavy
 */
bool Database::getUseLightProtein()
{
  return use_light_protein_;
}

/**
 *\returns the total number of proteins of the database
 */
unsigned int Database::getNumProteins()
{
  return proteins_->size();
}

/**
 *\returns the src FILE* of the database
 */
FILE* Database::getFile()
{
  return file_;
}

/**
 * sets the src FILE* of the database
 */
void Database::setFile(
  FILE* file ///< the src file to add -in
  )
{
  file_ = file;
}

/**
 * \returns the nth protein of the database
 * 
 */
Protein* Database::getProteinAtIdx(
  unsigned int protein_idx ///< The index of the protein to retrieve -in
  )
{
  //carp(CARP_DETAILED_DEBUG, "Getting db protein idx = %i, num proteins %i", 
  //     protein_idx, database->proteins.size());
  if( protein_idx >= proteins_->size()){
    carp(CARP_FATAL, 
         "Protein index %i out of bounds.  %i proteins in the database",
         protein_idx, proteins_->size());
  }

  return proteins_->at(protein_idx);
}

/**
 *\returns the protein designated by protein id of the database
 */
Protein* Database::getProteinByIdString(
  const char* protein_id ///< The id string for this protein -in
  ) {

  //TODO - Implement as a hashtable rather than a map to make 
  //this even faster if needed.
  Protein* protein = NULL;
  if (is_hashed_) {
    map<char*, Protein*, cmp_str>::iterator find_iter;
    find_iter = protein_map_->find((char*)protein_id);

    if (find_iter != protein_map_->end()) {
      protein = find_iter->second;
    }
  } else {
    //create the hashtable of protein ids
    for (unsigned int protein_idx = 0;
      protein_idx < proteins_->size();
      protein_idx++) {

      Protein* current_protein = proteins_->at(protein_idx);
      char* current_id = current_protein->getIdPointer();
      protein_map_->insert(make_pair(current_id, current_protein));

      if (strcmp(current_id, protein_id)==0) {
        protein = current_protein;
      }
        
    }
    is_hashed_ = true;
  }
  return protein;
}

/**
 * increase the pointer_count produced by this database.
 * \returns database pointer
 */
Database* Database::copyPtr(
  Database* database ///< the query database -in/out
  )
{
  if( database == NULL ){
    return NULL;
  }
  ++database->pointer_count_;
  return database;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

