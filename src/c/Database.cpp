/*************************************************************************//**
 * \file Database.cpp
 * \brief Object for representing a database of protein sequences.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
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
#include "sorter.h"
#include "ProteinIndex.h"
#include "ProteinIndexIterator.h"

#include "DatabaseProteinIterator.h"
#include "DatabasePeptideIterator.h"
#include "DatabaseSortedPeptideIterator.h"

#include <map>
#include <vector>

using namespace std;

/**
 * intializes a database object
 */
void Database::init(){
  filename_ = NULL;
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
  // fprintf(stderr, "Free: Allocation: %i\n", database->pointer_count);
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
  bool is_memmap ///< are we using a memory mapped binary fasta file? 
  ///< If so, all proteins are memory mapped -in
  )         
{
  carp(CARP_DEBUG, "Creating new database from '%s'", filename);
  init();
  setFilename(filename);
  is_memmap_ = is_memmap;
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

  free(filename_);
  
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
      if(munmap(data_address_, file_size_) != 0){
        carp(CARP_ERROR, "failed to unmap the memory of binary fasta file");
      }
    }
    // not memory mapped
    else{
      // close file handle
      carp(CARP_DEBUG, "Closing database filehandle");
      fclose(file_);
    }
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

  fprintf(file, "filename:%s\n", filename_);
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

  carp(CARP_DEBUG, "Parsing text fasta file '%s'", filename_);
  // check if already parsed
  if(is_parsed_){
    return true;
  }
  
  // open file and 
  file = fopen(filename_, "r");
  
  // check if succesfully opened file
  if(file == NULL){
    carp(CARP_ERROR, "Failed to open fasta file %s", filename_);
    return false;
  }
  
  // check if use light protein and parse thos light proteins fomr protein index
  if(use_light_protein_ && ProteinIndex::onDisk(filename_, false)){
    // let the user know that protein index file is being used
    carp(CARP_INFO, "using protein index file");

    // create a protein index iterator
    ProteinIndexIterator* protein_index_iterator =
      new ProteinIndexIterator(filename_);

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
            carp(CARP_ERROR, "failed to parse fasta file");
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
  if (stat(filename_, &file_info) == -1) {
    carp(CARP_ERROR,
         "Failed to retrieve information of binary fasta file: %s",
         filename_);
    return false;
  }
  
  // set size of the binary fasta file in database
  // this is used later to know how much to unmap
  file_size_ = file_info.st_size;
  
  // memory map the entire binary fasta file!
  data_address_ = mmap((caddr_t)0, file_info.st_size, PROT_READ, MAP_PRIVATE /*MAP_SHARED*/, file_d, 0);

  // check if memory mapping has succeeded
  if ((caddr_t)(data_address_) == (caddr_t)(-1)){
    carp(CARP_ERROR, "failed to use mmap function for binary fasta file: %s", filename_);
    return false;
  }
  
  return true;
}

/**
 * Assumes that there is a 1 at the very end after all the proteins in binary file
 *\return true successfully populates the proteins from memory mapped binary fasta file, else false
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
      carp(CARP_ERROR, "failed to parse fasta file");
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
  carp(CARP_DEBUG, "Parsing binary fasta file '%s'", filename_);
 
  // check if already parsed
  if(is_parsed_){
    return true;
  }
  
  // open file and 
  file_d = open(filename_, O_RDONLY);
  
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
  char* output_dir
  ){

  bool success = false;

  // from output_dir name and database filename, get binary_filename
  char* binary_filename = generate_name_path( filename_, ".fasta",
                                           "-binary-fasta", output_dir);


  carp(CARP_DEBUG, "Transforming text file '%s' to binary file '%s'",
       filename_, binary_filename);

  // create binary fasta
  success = create_binary_fasta_here(filename_, 
                                     binary_filename);

  if(! success ){
    carp(CARP_ERROR, 
         "Could not create binary fasta file '%s' from text fasta file '%s'", 
         binary_filename, filename_);
    return false;
  }
  // change database filename to new binary fasta
  char* binary_filename_no_path = parse_filename(binary_filename);
  setFilename(binary_filename);

  // set is_memmap to true
  setMemmap(true);

  // parse the binary fasta
  success = parse();

  free(binary_filename);
  free(binary_filename_no_path);
  return success;
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
  return my_copy_string(filename_);
}


/**
 *\returns the pointer to the filename of the database
 * user must not free or change the filename
 */
char* Database::getFilenamePointer()
{
  return filename_;
}

/**
 * sets the filename of the database
 * protein->sequence must been initiailized
 */
void Database::setFilename(
  const char* filename ///< the filename to add -in
  )
{
  free(filename_);
  // MEMLEAK below
  filename_ = my_copy_string(filename);
}

/**
 *\returns true|false whether the database has been parsed?
 */
bool Database::getIsParsed()
{
  return is_parsed_;
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
    map<char*, Protein*>::iterator find_iter;
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




/***********************************************
 * Iterators
 ***********************************************/




/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_peptide_iterator object.
 */
void void_free_database_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  )
{

  delete (DatabasePeptideIterator*)database_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns true if there are additional peptides, false if not.
 */
bool void_database_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return ((DatabasePeptideIterator*)database_peptide_iterator)->hasNext();
}

/**
 * \returns The next peptide in the database.
 */
Peptide* void_database_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{

  return ((DatabasePeptideIterator*)database_peptide_iterator)->next();
}

/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 * Frees an allocated database_sorted_peptide_iterator object.
 */
void void_free_database_sorted_peptide_iterator(
  void* database_peptide_iterator ///< the iterator to free -in
  )
{
  delete (DatabaseSortedPeptideIterator*)database_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns true if there are additional peptides to iterate over, false if not.
 */
bool void_database_sorted_peptide_iterator_has_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return ((DatabaseSortedPeptideIterator*)database_peptide_iterator)->hasNext();
}

/**
 * returns each peptide in sorted order
 * \returns The next peptide in the database.
 */
Peptide* void_database_sorted_peptide_iterator_next(
  void* database_peptide_iterator ///< the iterator of interest -in
  )
{
  return ((DatabaseSortedPeptideIterator*)database_peptide_iterator)->next();
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

