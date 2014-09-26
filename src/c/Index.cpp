/************************************************************************//**
 * \file index.cpp
 * \brief Object for representing an index of a database
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#endif
#include <time.h>
#include <errno.h>
#include <iostream>
#include <sstream>

#include "utils.h"
#include "crux-utils.h"

#include "Peptide.h"
#include "Protein.h"
#include "Index.h"
#include "IndexPeptideIterator.h"
#include "IndexFilteredPeptideIterator.h"
#include "carp.h"
#include "sorter.h"
#include "objects.h"
#include "PeptideConstraint.h"
#include "Database.h"
#include "DatabasePeptideIterator.h"
#include "ProteinIndex.h"
#include "parameter.h"
#include "WinCrux.h"

using namespace std;
using namespace Crux;

// maximum proteins the index can handle
static const int MAX_PROTEIN = 30000;
static const int NUM_CHECK_LINES = 8;
static const int MAX_PROTEIN_IN_BIN = 2500;
static const int MAX_FILE_SIZE_TO_USE_LIGHT_PROTEIN = 500000000;
const char* Index::index_file_prefix = "crux_index_";
const char* Index::decoy_index_file_prefix = "crux_decoy_index_";


#ifdef DARWIN
// OS X doesn't support fcloseall()
// Consider tracking open file handles
// rather then simply relying on process
// termination to close files.
#define fcloseall()
#endif

/* 
 * How does the create_index (the routine) work?
 *
 * - Creates a database object from a the fasta file member variable.
 *
 * - Depending on the size of the database, determines how many passes it
 *   will need to make over the database peptides (first implementation
 *   will only make one pass, later implementation can make multiple passes
 *   with multiple iterators)
 *
 * - From the directory and mass resolution member variables, 
 *   creates the list of filenames necessary for storing the peptides 
 *   
 * - Creates a database_peptide_iterator object from the database 
 *   and peptide constraint objects
 *
 * - Then starts iterating through peptides 
 *
 *    - From the peptide mass determines which file to place it in, and
 *      then calls the serialize_peptide method
 *
 *      - The serialize peptide method is an object method specifically
 *        designed for this purpose, and writes a peptide in a format from
 *        which it can be reconstructed
 *
 *        - Serialize peptide needs to serialize the peptide_src objects. To
 *        do that, it needs the idx of the peptide_src's protein objects in
 *        the database. This is retrieved from the protein idx member
 *        variable (i.e. this protein is the 0th protein, 1st protein etc.)
 *        which is set at database creation. Note, that this won't have to 
 *        change when we move to light proteins.
 *
 * LATER At some point we will index each of the peptide files too (in
 * a TOC), so that we can retrieve them rapidly later. 
 *
 * LATER We implement light proteins and possibly an 
 * LATER create_index for protein locations in the database allowing rapid
 * creation of light protein objects.
 */

/* 
 * How does generate_peptides (the executable) with --from-index work?
 *
 * - Given a fasta file, looks for the corresponding index on disk. Fail if it
 *   can't find the corresponding index files.
 *
 * - Instantiates an index object from the fasta file with new_index.
 *
 * - Attempts to parse_index the index.
 *
 * - Instantiates an index_peptide_iterator from the index, according to
 *   constraints passed on the command line to generate_peptides.
 *
 * - Then starts iterating through peptides, which are being loaded from
 *   disk, and outputs them as before
 * 
 * LATER We implement light proteins.
 * LATER use an index for protein offsets in the database allowing rapid
 * creation of light protein objects.
 *
 * LATER Develop a conversion from light to heavy and heavy to light protein
 * objects to avoid excessive memory use.
 */

/**
 * \struct bin_peptide_iterator
 * \brief An iterator to iterate over the peptides in a bin (one file
 * handle) 
 */
struct bin_peptide_iterator{
  Index* index; ///< The index object which we are iterating over
  Database* database; ///< either target or decoy of index
  FILE* index_file; ///< The current file stream that we are reading from
  bool has_next; ///< Is there another peptide?
  Peptide* peptide; ///< the next peptide to return
};    


/**
 * \struct bin_sorted_peptide_iterator
 * \brief Object to iterate over the peptides within a bin, in an
 * sort in mass
 */
struct bin_sorted_peptide_iterator {
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator; 
  ///< the peptide iterator that sorts the peptides form the bin
};

/************
 * Public functions for Index
 ************/


/**
 * initializes an Index object.
 */
void Index::init() {
  num_pointers_ = 1;
  database_ = NULL;
  decoy_database_ = NULL;
  index_map_ = NULL;
  decoy_index_map_ = NULL;
  directory_ = NULL;
  disk_constraint_ = NULL;
  search_constraint_ = NULL;
  on_disk_ = false;
  mass_range_ = 0.0;
  is_unique_ = false;
  decoys_ = NO_DECOYS;
  static_mods_ = new FLOAT_T[Alphabet::numAminoAcids];
  for(int aa_idx = 0; aa_idx < Alphabet::numAminoAcids; aa_idx++){
    static_mods_[aa_idx] = 0;
  }
}  

/**
 * \returns An (empty) index object.
 */
Index::Index() {
  init();
}

/**
 * Helper function for scandir to find files with names ending in
 * Database::binary_suffix.
 * \returns 1 if filename has ending, else 0
 */
#if ((defined DARWIN && defined PRE_MOUNTAIN_LION) || defined WIN32)
int is_binary_fasta_name(struct dirent *entry){
#else
int is_binary_fasta_name(const struct dirent *entry){
#endif


  const char* filename = entry->d_name; //w/o const gcc warning
  const char* suffix = Database::binary_suffix.c_str();

  int name_length = strlen(filename);
  int suffix_length = strlen(suffix);

  if( suffix_length > name_length){
    return 0;
  }

  // point to the last bit of filename
  filename += name_length - suffix_length;
  int matches = strncmp(filename, suffix, suffix_length); // 0 if matches

  if( matches == 0 ) {
    return 1;
  }//else

  return 0;
}

/**
 * Helper function for scandir to find files with names ending in
 * Database::decoy_binary_suffix.
 * \returns 1 if filename has ending, else 0. 
 */
#if ((defined DARWIN && defined PRE_MOUNTAIN_LION) || defined WIN32)
int is_decoy_binary_fasta_name(struct dirent *entry){
#else
int is_decoy_binary_fasta_name(const struct dirent *entry){
#endif


  const char* filename = entry->d_name; //w/o const gcc warning
  const char* suffix = Database::decoy_binary_suffix.c_str();

  int name_length = strlen(filename);
  int suffix_length = strlen(suffix);

  if( suffix_length > name_length){
    return 0;
  }

  // point to the last bit of filename
  filename += name_length - suffix_length;
  int matches = strncmp(filename, suffix, suffix_length); // 0 if matches

  if( matches == 0 ) {
    return 1;
  }//else

  return 0;
}

/**
 * \brief Looks in given directory for a file ending in
 * Database::binary_suffix and returns a heap-allocated string of the full
 * name including the index directory.
 *
 * Exits with error if index->directory does not exist, no file
 * exists, or more than one file exists. 
 * \returns A string with the name of the existing binary fasta file
 * for this index. 
 */
char* Index::getBinaryFastaName(
 const char* index_name ///< directory containing index
){
  struct dirent** namelist = NULL;
  int num_files = scandir(index_name, &namelist, is_binary_fasta_name,
                          alphasort);

  if( num_files < 1 ){
    carp(CARP_FATAL, "Binary fasta file missing from index '%s'", 
         index_name);
  }
  if( num_files > 1 ){
    carp(CARP_FATAL, "More than one binary fasta file found in index '%s'", 
         index_name);
  }

  carp(CARP_DEBUG, "Found '%s' in index '%s'", namelist[0]->d_name, 
       index_name);

  char* fasta_name = my_copy_string(namelist[0]->d_name);
  char* fasta_name_path = get_full_filename(index_name, fasta_name);

  std::free(namelist[0]);
  std::free(namelist);
  std::free(fasta_name);

  return fasta_name_path;

}

/**
 * \brief Looks in given directory for a file ending in
 * Database::decoy_binary_suffix and returns a heap-allocated string of the
 * full name including the index directory.
 *
 * Exits with error if index_name does not exist or more than one file exists.
 * \returns A string with the name of the existing binary fasta file
 * for this index or NULL if no file found.
 */
char* Index::getDecoyBinaryFastaName(
  const char* index_name ///< directory containing index
){
  struct dirent** namelist = NULL;
  int num_files = scandir(index_name, &namelist, is_decoy_binary_fasta_name,
                          alphasort);

  if( num_files < 1 ){
    return NULL;
  }
  if( num_files > 1 ){
    carp(CARP_FATAL, "More than one decoy binary fasta file found in '%s'", 
         index_name);
  }

  carp(CARP_DEBUG, "Found '%s' in index '%s'", namelist[0]->d_name, 
       index_name);

  char* fasta_name = my_copy_string(namelist[0]->d_name);
  char* fasta_name_path = get_full_filename(index_name, fasta_name);

  std::free(namelist[0]);
  std::free(namelist);
  std::free(fasta_name);

  return fasta_name_path;

}
/**
 * \brief Initializes a new Index object by setting its member
 * variables to the values in the index on disk.
 *
 * Assumes that the name of the index (directory) has been set by
 * caller.  Checks that the directory exists.  Reads the mapfile
 * therein and gets the constraint values from the header.  Returns
 * true if all of the above are true, else prints an error and returns
 * false. 
 *
 * \returns true if fields were successfully set, false if index on
 * disk could not be used to initialize the object.
 */
bool Index::setFieldsFromDisk() {

  assert(directory_ != NULL );
  DIR* check_dir = opendir(directory_);
  
  if(check_dir == NULL){
    carp(CARP_ERROR, "Failed to access index directory %s", directory_);
    return false;
  }else{
    closedir(check_dir);
  }

  // create constraint
  disk_constraint_ = new PeptideConstraint();
  search_constraint_ = NULL;

  // open map file
  char* map_filename = get_full_filename(directory_, "crux_index_map");
  carp(CARP_DEBUG, "Opening index map file '%s'", map_filename);
  FILE* map_file = fopen(map_filename, "r");
  if(map_file == NULL){
    carp(CARP_ERROR, "Could not open index map file %s", map_filename);
    return false;
  }

  // read header lines (beginning with #)
  char* line = NULL;
  size_t buf_length = 0;
  int line_length = getline(&line, &buf_length, map_file);
  while( line_length != -1 && line[0] == '#' ){
    setFieldFromMap(line);
    line_length = getline(&line, &buf_length, map_file);
  }
  fclose(map_file);
  std::free(map_filename);
  myfree(line);

  // check that the static mods in the index are the same as the parameters
  for(int aa_idx = 0; aa_idx < Alphabet::numAminoAcids; aa_idx++){
    const char* amino_acid = Alphabet::aminoAcids[aa_idx];
    FLOAT_T param_mod = get_double_parameter(amino_acid);
    if( static_mods_[aa_idx] != param_mod ){
      carp(CARP_FATAL, "Index cannot be searched with static mod "
           "%s=%f.  The index's static mod for %s is %f.",
           amino_acid, param_mod, amino_acid, static_mods_[aa_idx]);
    }
  }

  on_disk_ = true;
  return true;
}

/**
 * \brief Private function to take a header line from index_map and
 * set the appropriate value in the index.
 */
void Index::setFieldFromMap(
  char* line ///< line read from map file
){
  // parse the line
  char sharp[2] = "";
  char trait_name[64] = "";
  FLOAT_T value = 0;
  #ifdef USE_DOUBLES
  sscanf(line, "%s %s %lf", sharp, trait_name, &value);
  #else
  sscanf(line, "%s %s %f", sharp, trait_name, &value);
  #endif
  carp(CARP_DEBUG, "Index map header found %s value %.2f", trait_name, value);

  if(strcmp("min_mass:", trait_name) == 0){
    disk_constraint_->setMinMass(value);
  }
  else if(strcmp("max_mass:", trait_name) == 0){
    disk_constraint_->setMaxMass(value);
  }
  else if(strcmp("min_length:", trait_name) == 0){
    disk_constraint_->setMinLength((int)value);
  }
  else if(strcmp("max_length:", trait_name) == 0){
    disk_constraint_->setMaxLength((int)value);
  }
  else if(strcmp("peptide_type:", trait_name) == 0){
    //    set_peptide_constraint_peptide_type(index->disk_constraint, value);
    carp(CARP_FATAL, "This index uses the obsolete 'peptide_type'. "
         "Rebuild with current version of crux to use 'enzyme_type'.");
  }
  else if(strcmp("enzyme_type:", trait_name) == 0){
    disk_constraint_->setEnzyme((ENZYME_T) ((int) value));
  }
  else if(strcmp("digest_type:", trait_name) == 0){
    disk_constraint_->setDigest(((DIGEST_T) (int) value));
  }
  else if(strcmp("missed_cleavages:", trait_name) == 0){
    disk_constraint_->setNumMisCleavage((int)value);
  }
  else if(strcmp("mass_type:", trait_name) == 0){
    disk_constraint_->setMassType((MASS_TYPE_T)((int) value));
  }
  else if(strcmp("unique_peptides:", trait_name) == 0){
    is_unique_ = (bool)value;
  }
  else if(strcmp("target_mass_range_for_index_file:", trait_name) == 0){
    mass_range_ = value;
  }
  else if(strcmp("decoys:",trait_name) == 0 ){
    decoys_ = (DECOY_TYPE_T)((int) value);
  }
  else if(strncmp("static_mod_", trait_name, strlen("static_mod_")) == 0){
    const char* aa = trait_name + strlen("static_mod_");
    int aa_idx = Alphabet::aminoToInt(*aa);
    static_mods_[aa_idx] = value;
  }
  else if(strcmp("CRUX", trait_name) == 0){
    return;
  }
  else if(strcmp("time", trait_name) == 0){
    return;
  }
  else{
    carp(CARP_FATAL, "Unknown index map entry: %s (%s, %d)", 
         line, trait_name, value);
  }
}

/**
 * \brief Private function to confirm that the attributes of the index
 * will work for this run.
 * Requires that the range from min to max length and mass include the
 * range requested in parameter.c.  The index must have at least as
 * many missed cleavages as requested.  The mass type must be the
 * same. Requires index and request be both unique or both not 
 * unique.  The peptide cleavage type of the index must be no more
 * restrictive than requested (where TRYPTIC is most restrictive and
 * ANY_TRYPTIC is least).
 * \returns true if all constraints will work with those in parameter.c.
 */
bool Index::checkConstraints(){
  double min_mass = disk_constraint_->getMinMass();
  double max_mass = disk_constraint_->getMaxMass();
  double min_len = disk_constraint_->getMinLength();
  double max_len = disk_constraint_->getMaxLength();
  int missed_cleavages = 
    disk_constraint_->getNumMisCleavage();
  ENZYME_T enzyme = disk_constraint_->getEnzyme();
  DIGEST_T digestion = disk_constraint_->getDigest();
  MASS_TYPE_T mass_type = disk_constraint_->getMassType();

  bool success = true;
  const char* param;
  if(min_mass > get_double_parameter("min-mass")){
    success = false;
    param = "min-mass";
  }else if(max_mass < get_double_parameter("max-mass")){
    success = false;
    param = "max-mass";
  }else if(min_len > get_int_parameter("min-length")){
    success = false;
    param = "min-length";
  }else if(max_len < get_int_parameter("max-length")){
    success = false;
    param = "max-length";
  }else if(missed_cleavages < get_int_parameter("missed-cleavages")){
    success = false;
    param = "missed-cleavages";
  }else if(mass_type != get_mass_type_parameter("isotopic-mass")){
    success = false;
    param = "isotopic-mass";
  }else if(enzyme != get_enzyme_type_parameter("enzyme")){
    printf("constraint e: %i, param e: %i.", (int)enzyme, (int)get_enzyme_type_parameter("enzyme"));
    success = false;
    param = "enzyme";
  }else if(digestion < get_digest_type_parameter("digestion")){
    // full < partial < non-specific
    success = false;
    param = "digestion";
  }

  if(success == false){
    carp(CARP_ERROR, "Index does not support the given search parameters, " \
         "'%s' is not compatable", param);
  }
  return success;
}

/**
 * \brief The initialization step in creating an index after
 * allocation and before parsing database.  Private function called
 * only by new_index().  No assumptions about fasta
 * or index filenames or locations relative to eachother.  Does not
 * initialize the database to be associated with this index.
 *
 * \returns void
 */
void Index::setFields(
  const char* output_dir,      ///< The name of the new index
  PeptideConstraint* constraint,  
  ///< Constraint which these peptides satisfy -in
  FLOAT_T mass_range,  
  ///< the range of mass that each index file should be partitioned into -in
  bool is_unique, ///< only unique peptides? -in
  DECOY_TYPE_T decoys ///< the type of decoys to store
  )
{
  carp(CARP_DEBUG, "Setting index fields");
  DIR* check_dir = NULL;
  
  if((check_dir = opendir(output_dir)) != NULL){
    on_disk_ = true;
    closedir(check_dir);
  }
  else{
    on_disk_ = false;
  }

  carp(CARP_DEBUG, "Index on disk is '%i' for dir '%s'", 
       (int)on_disk_, output_dir);

  // set each field
  setDirectory(output_dir);
  disk_constraint_ = constraint;
  search_constraint_ = NULL;
  //index->on_disk = false; // this breaks overwrite of create index
  mass_range_ = mass_range;  
  is_unique_ = is_unique;
  decoys_ = decoys;
}


/**
 * \brief Create a new index for the given fasta file of peptides
 * conforming to the given constraints.
 *
 * Allocates memory for the index.  Allocates memory for the database
 * (object for reading the fasta file) and associates it with the
 * index. Does not read in peptides from the fasta file.
 * Assumes that the fasta file is always in the directory where the
 * crux_index_file directory is located.
 * USE this constructor when creating index, does not parse database 
 * database is later transformed into memory mapped database and parse
 * in create index.  For creating peptides from an index on disk, use
 * new_index_from_disk routine 
 * \returns A new index object.
 */
Index::Index(
  const char* fasta_filename, ///< The fasta file
  const char* output_dir,     ///< The name of the new index
  PeptideConstraint* constraint,  
    ///< Constraint which these peptides will satisfy
  FLOAT_T mass_range,  
    ///< the range of mass that each index file should be partitioned into
  DECOY_TYPE_T decoys ///< what kind of decoys to produce, if any
  ) {

  carp(CARP_DETAILED_DEBUG, "Creating new index to be named %s", output_dir);
  init();
  Database* database = NULL;
  
  // Initially, create a database that does not use memory mapping (false)
  // Once binary fasta file has been creaated this will change to a
  // memory mapped database
  database = new Database(fasta_filename, false);

  // set database, has not been parsed
  database_ = database;
  if( decoys == NO_DECOYS ){
    decoy_database_ = NULL;
  } else {
    decoy_database_ = new Database(fasta_filename, false, decoys);
  }
  bool is_unique = true;
  setFields(output_dir, constraint, mass_range, is_unique, decoys);
}         

/**
 * \brief Create an index object to represent an index already on disk.
 * Used by generate-peptides and search-for-matches.  If the proper
 * files are not found in the directory, warns and returns NULL.
 *
 * \returns A new index object ready for search.
 */
Index::Index(
  const char* index_name  ///< The directory containing the index
  ) {

  // find the database file, open it and point to it
  // check that it works with global constraints

  // allocate index and name it
  init();
  setDirectory(index_name);
  
  // this should find the map file and get constraint values
  if(!setFieldsFromDisk()){
    carp(CARP_FATAL, "Could not create index %s from disk", index_name);
  }

  // check that the index constraints are OK for this run
  if(!checkConstraints()){
    carp(CARP_FATAL, "Index cannot produce the requested peptides. " \
     "Change the search parameters to match the index or create new index.");
  }

  // get binary fasta file name with path to crux directory 
  char* binary_fasta = getBinaryFastaName(index_name);
  
  // check if input file exist
  if(access(binary_fasta, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist for crux index.", 
        binary_fasta);
  }
  
  // now create a database, using binary fasta file
  database_ = new Database(binary_fasta, true);
  
  if(!database_->parse()){
    carp(CARP_FATAL, "Failed to parse database, cannot create new index");
  }

  // free string
  std::free(binary_fasta);

  // repeat for decoy database
  if( decoys_ != NO_DECOYS ){
    char* decoy_binary_fasta = getDecoyBinaryFastaName(index_name);

    // check if input file exist
    if(access(decoy_binary_fasta, F_OK)){
      carp(CARP_FATAL, "The file \"%s\" does not exist for crux index.", 
           decoy_binary_fasta);
    }
    decoy_database_ = new Database(decoy_binary_fasta, true, decoys_);
    
    if(!decoy_database_->parse()){
      carp(CARP_FATAL, "Failed to parse decoy database, can't create index");
    }

    // free string
    std::free(decoy_binary_fasta);
  }
}

/**
 * \returns The number of proteins in the index
 */
int Index::getNumProteins(){
  // will be the same for decoy database, if present
  return database_->getNumProteins();
}

/**
 * Increase the number of objects pointing to this.
 */
Index* Index::copyPtr() {

  num_pointers_++;
  return this;
}

/**
 * Frees an allocated index object, but only if pointer count is equal to 1.
 */
void Index::free(
  Index* index
  )
{
  if(index == NULL ){
    return;
  }

  if (index->num_pointers_ > 1){
    index->num_pointers_--;
  } else {
    carp(CARP_DEBUG, "Freeing index");
    delete index;
  }
}

/**
 * Destructor that frees the indexes database(s) and constraints.
 */
Index::~Index() {

  if (database_ != NULL){
    carp(CARP_DEBUG, "Freeing index database");
    Database::freeDatabase(database_);
  }

  if (decoy_database_ != NULL){
    carp(CARP_DEBUG, "Freeing index decoy database");
    Database::freeDatabase(decoy_database_);
  }

  if (disk_constraint_ != NULL){
    carp(CARP_DEBUG, "Freeing index disk constraint");
    PeptideConstraint::free(disk_constraint_);
  }

  if (search_constraint_ != NULL){
      carp(CARP_DEBUG, "Freeing index search constraint");
      PeptideConstraint::free(search_constraint_);
    }

  std::free(directory_);
  delete [] static_mods_;
}

/**
 * Write to the file stream various information of the
 * index files created.
 */
bool Index::writeHeader(
  FILE* file ///< out put stream for crux_index_map -in
  )
{
  time_t hold_time;
  hold_time = time(0);
  PeptideConstraint* constraint = disk_constraint_;
  
  fprintf(file, "#\tmin_mass: %.2f\n", constraint->getMinMass());
  fprintf(file, "#\tmax_mass: %.2f\n", constraint->getMaxMass());
  fprintf(file, "#\tmin_length: %d\n", constraint->getMinLength());
  fprintf(file, "#\tmax_length: %d\n", constraint->getMaxLength());
  fprintf(file, "#\tenzyme_type: %d\n", constraint->getEnzyme());
  fprintf(file, "#\tdigest_type: %d\n", constraint->getDigest());
  fprintf(file, "#\tmissed_cleavages: %d\n", constraint->getNumMisCleavage());
  fprintf(file, "#\tmass_type: %d\n", constraint->getMassType());
  fprintf(file, "#\tunique_peptides: %d\n", getIsUnique());
  fprintf(file, "#\tdecoys: %d\n", decoys_);
  
  fprintf(file, "#\tCRUX index directory: %s\n", directory_);
  fprintf(file, "#\ttime created: %s",  ctime(&hold_time)); 
  fprintf(file, "#\ttarget_mass_range_for_index_file: %.2f\n", mass_range_);
  
  for(int aa_idx = 0; aa_idx < Alphabet::numAminoAcids; aa_idx++){
    const char* amino_acid = Alphabet::aminoAcids[aa_idx];
    FLOAT_T mod =  get_double_parameter(amino_acid);
    if( mod != 0 ){
      fprintf(file, "#\tstatic_mod_%s: %+.6f\n", amino_acid, mod);
    }
  }
  return true;
}

/**
 * write to the file stream various information of the
 * index files created in human readable format
 *\returns true if successfully creates README file, else false
 */
bool Index::writeReadmeFile(
  FILE* file ///< out put stream for README file -in
  )
{
  time_t hold_time;
  hold_time = time(0);
  PeptideConstraint* constraint = disk_constraint_;
  char* fasta_file = database_->getFilename();
  char* fasta_file_no_path = parse_filename(fasta_file);
  char* decoy_fasta_file_no_path = NULL;
  if( decoys_ != NO_DECOYS ){
    char* decoy_fasta_file = decoy_database_->getFilename();
    decoy_fasta_file_no_path = parse_filename(decoy_fasta_file);
    std::free(decoy_fasta_file);
  }
  
  fprintf(file, "#\ttime created: %s",  ctime(&hold_time)); 
  fprintf(file, "#\tfasta file: %s\n",  fasta_file_no_path); 
  if( decoys_ != NO_DECOYS ){
    fprintf(file, "#\tdecoy fasta file: %s\n",  decoy_fasta_file_no_path); 
  }
  fprintf(file, "#\tmin_mass: %.2f\n",
          constraint->getMinMass());
  fprintf(file, "#\tmax_mass: %.2f\n",
          constraint->getMaxMass());
  fprintf(file, "#\tmin_length: %d\n",
          constraint->getMinLength());
  fprintf(file, "#\tmax_length: %d\n",
          constraint->getMaxLength());

  char* value_str = 
    enzyme_type_to_string(constraint->getEnzyme());
  fprintf(file, "#\tenzyme_type: %s\n", value_str);
  std::free(value_str);
  value_str = 
    digest_type_to_string(constraint->getDigest());
  std::free(value_str);

  fprintf(file, "#\tmissed_cleavage: %s\n", 
        (constraint->getNumMisCleavage()? "true":"false"));
  fprintf(file, "#\tmass_type: %s\n", 
          (constraint->getMassType()==AVERAGE? 
           "average":"mono"));
  fprintf(file, "#\tunique peptides: %s\n", 
          (getIsUnique()? "unique":"redundant"));
  fprintf(file, "#\tCRUX index directory: %s\n", directory_);
  fprintf(file, "#\ttarget mass range for index file: %.2f\n", 
          mass_range_);
  
  std::free(fasta_file);
  std::free(fasta_file_no_path);
  std::free(decoy_fasta_file_no_path);

  return true;
}


/**
 * Combines the bin index and the file prefix into the name of the
 * index file for that bin.
 *\returns the filename for the given index
 */
char* get_crux_filename(
  long bin_idx,  ///< the bin_indx name you want 
  const char* file_prefix  ///< the first part of the filename
  )
{
  ostringstream filename;
  filename << file_prefix;
  filename << (bin_idx + 1);

  return my_copy_string(filename.str().c_str());
}

/**
 * Calculate the total number of bins( file handlers) that will be needed
 * \returns the total number of bins needed
 */
long Index::getNumBinsNeeded(
  int* mass_limits  ///< an array that holds the min/max mass limit -in
  )
{
  int min_length = disk_constraint_->getMinLength();
  int max_length = disk_constraint_->getMaxLength();
  FLOAT_T min_mass = disk_constraint_->getMinMass();
  FLOAT_T max_mass = disk_constraint_->getMaxMass();
  FLOAT_T min_mass_limit = min_mass;
  FLOAT_T max_mass_limit = max_mass;
  long num_bins = 0;

  // reset minimum mass limit
  if(min_length * 57 + MASS_H2O_MONO > min_mass){
    min_mass_limit = min_length * 57 + MASS_H2O_MONO; 
  }
  
  // reset maximum mass limit
  if(max_length * 187 + MASS_H2O_AVERAGE < max_mass){
    max_mass_limit = max_length * 187 + MASS_H2O_AVERAGE;
  }

  // set mass limit info array
  min_mass_limit = (int)min_mass_limit;
  max_mass_limit = (int)max_mass_limit + 1;

  (mass_limits)[0] = (int)min_mass_limit;
  (mass_limits)[1] = (int)max_mass_limit;

  num_bins = (int)(((max_mass_limit - min_mass_limit) / mass_range_) + 1);  // check..

  if (num_bins > MAX_INDEX_FILES){
    num_bins = MAX_INDEX_FILES;
  }
  return num_bins;
}                         

/**
 * Opens one index file, storing it in the array.  File name is the
 * concatenation of file_prefix and bin_index+1.
 * \returns True, if successfully opened file, else false.
 */
bool generate_one_file_handle(
  FILE** file_array,  ///< the file handler array 
  long bin_index, ///< the bin index for which to create a file handle
  const char* file_prefix
  )
{
  // create the file handler needed for create index
  char* filename = get_crux_filename(bin_index, file_prefix);
  FILE* file = fopen(filename, "w+" );
  free(filename);
  if(file == NULL){
    carp(CARP_WARNING, "cannot open all file handlers needed");
    return false;
  }
  file_array[bin_index] = file;
  return true;
}

/**
 * given a mass, it finds the correct bin for that mass and returns
 * the file handler to that bin
 *\returns the FILE handler that the mass should be collected
 */
FILE* get_bin_file(
  int mass, ///< the query mass -in
  int low_mass,///< the mass limit array -in
  FILE*** file_array ///< file handler array -in
  )
{
  FILE* file = NULL;
  long bin_idx = (mass - low_mass)/100 + 1;
  file = (*file_array)[bin_idx];
  return file;
}

/**
 * \brief For one file bin, reparses the peptides, sorts them, and
 * reprints them to the crux_index file 
 *
 * \returns the sorted bin
 */
FILE* Index::sortBin(
  FILE* file, ///< the working file handle to the bin -in
  long bin_idx, ///< bin index in the file array -in
  unsigned int peptide_count, ///< the total peptide count in the bin -in
  FILE* text_file, ///< optional file to print sequences to -in
  const char* file_prefix, ///< beginning of index file names -in
  Database* database ///< target or decoy database to use -in
  )
{
  char* filename = NULL;

  // check if file is empty
  if(peptide_count == 0){
    return file;
  }

  BIN_SORTED_PEPTIDE_ITERATOR_T* peptide_iterator =
    new_bin_sorted_peptide_iterator(this, file, peptide_count, database);
  Peptide* working_peptide = NULL;
  
  // get the filename for this file bin
  filename = get_crux_filename(bin_idx, file_prefix);
  // close unsorted bin
  fclose(file);
  // create new bin which will be sorted 
  file = fopen(filename, "w" );
  
  // serialize all peptides in sorted order
  while(bin_sorted_peptide_iterator_has_next(peptide_iterator)){
    working_peptide = bin_sorted_peptide_iterator_next(peptide_iterator);
    working_peptide->serialize(file, text_file);
    delete working_peptide;
  }
  
  std::free(filename);
  free_bin_sorted_peptide_iterator(peptide_iterator);

  return file;
}

/**
 * Stores the peptide in the correct bin.  If bin exceeds
 * MAX_PROTEIN_IN_BIN, then serialize all peptides in the bin.
 *
 * \returns true if successful in storing the peptide or serializing
 * peptides, else false .
 */
static bool dump_peptide(
  FILE** file_array,   ///< the working file handler array to the bins -in/out
  long int file_idx,   ///< the index of the file the peptide belongs to -in
  Peptide* working_peptide, ///< the peptide to be stored -in
  Peptide** peptide_array, ///< hold peptides before they're serialized -out
  int* bin_count       ///< an array of the counts of peptides in each bin -in
  )
{  
  int peptide_idx = 0;
  FILE* file = NULL;
  int current_count;
  
  // if the peptide count is over the limit
  if((current_count = bin_count[file_idx]) >= MAX_PROTEIN_IN_BIN){
    file = file_array[file_idx];
    // print out all peptides
    for(peptide_idx = 0; peptide_idx < current_count; ++peptide_idx){
      peptide_array[peptide_idx]->serialize(file, NULL);
      delete peptide_array[peptide_idx];
      // FIXME this should not be allowed
    }
    working_peptide->serialize(file, NULL);
    delete working_peptide;
    // FIXME this should not be allowed
    bin_count[file_idx] = 0;
  }
  // if the peptide count is bellow the limit
  else{
    // store peptide in peptide array , these peptides will be printed later togehter
    peptide_array[(bin_count[file_idx])] = working_peptide;
    ++bin_count[file_idx];
  }
  return true;
}

/**
 * /brief Serializes all peptides in the peptide array after parsing
 * the database.  Used in creating index after parsing database and
 * before sorting peptides in each bin.
 *
 * As a side effect, if the last parameter is non-NULL, print the
 * peptides in ASCII format to the specified file.  Frees both the
 * peptide_array and bin_count once all serialized.
 *
 * \returns true, if successful in serializing all peptides, else false.
 */
static bool dump_peptide_all(
  FILE** file_array,   ///< the working file handle array to the bins -out
  Peptide*** peptide_array, ///< the array of pre-serialized peptides -in
  int* bin_count,      ///< the count array of peptides in each bin -in
  int num_bins         ///< the total number of bins -in
  )
{  
  int peptide_idx = 0;
  FILE* file = NULL;
  Peptide** working_array = NULL;
  int bin_idx = 0;
  int file_idx = 0;

  // print out all remaining peptides in the file_array
  for(file_idx = 0; file_idx < num_bins; ++file_idx){
    carp(CARP_DETAILED_DEBUG, "Serializing bin %d", file_idx);
    // no peptides in this bin
    if((file = file_array[file_idx]) == NULL){
      free(peptide_array[file_idx]);
      continue;
    }
    working_array = peptide_array[file_idx];
    // BF: could we use int bin_count = bin_count[file_idx]
    bin_idx = bin_count[file_idx];
    // print out all peptides in this specific bin
    // BF: could we do for(pep_idx=0; pep_idx<bin_count_total; pep_idx++)
    while(bin_idx > 0){
      working_array[peptide_idx]->serialize(file, NULL);
      delete working_array[peptide_idx];
      // FIXME this should not be allowed
      --bin_idx;
      ++peptide_idx;
    }
    peptide_idx = 0;
    free(peptide_array[file_idx]);
  }
  // free both peptide array and bin_count array , bye bye
  free(bin_count);
  free(peptide_array);
  return true;
}

/**
 * \brief The main index method. Does all the heavy lifting, creating
 * files serializing peptides, etc. 
 *
 * When this is called, the database should be allocated and
 * initialized, but not parsed.  The peptide constraint, mass range,
 * and unique fields have been set.  If --overwrite is true, the output
 * dir may exist and still contain an index.
 *
 * The index directory itself should have  
 * a standard suffix (e.g. cruxidx), so that a given fasta file will have
 * an obvious index location.
 *
 * assumes that the current directory is the crux directory where the 
 * fasta file is located
 *
 * Note: creates an index .info file as to map masses to files and to store 
 * all information that was used to create this index.
 * \returns true if success. false if failure.
 */
bool Index::create(
  bool create_text_file ///< Should an ASCII text file be create? -in
  )
{
  carp(CARP_DEBUG, "Creating index");

#ifdef _MSC_VER
  // Windows defaults to a fairly small number
  // for the max open files. Increase this to the
  // maximum allowed by the Microsoft CRT.
  _setmaxstdio(WIN_MAX_OPEN_FILES);
#endif
  
  if(on_disk_){
    if (get_boolean_parameter("overwrite")) {
      carp(CARP_DEBUG, "Will be replacing existing index");
    } else { // this should have already been checked, but...
      carp(CARP_FATAL, "Index '%s' already exists.  " \
      "Use --overwrite T to replace", directory_);
    }
  }
  
  // Create index directory if it doesn't already exist.
  if (!on_disk_) {
    if(mkdir(directory_, S_IRWXU+S_IRWXG+S_IROTH+S_IXOTH) != 0){
      carp(CARP_WARNING, "Cannot create index directory %s", directory_);
      return false;
    }
  }

  // instantiate database(s)
  if(!database_->transformTextToMemmap(directory_, false)){//binary not tmp
    carp(CARP_FATAL, "Failed to create binary database from text fasta file");
  }

  if( decoy_database_ 
      && ! decoy_database_->transformTextToMemmap(directory_, false) ){
    carp(CARP_FATAL, "Failed to create binary decoy database from fasta file");
  }

  // move into index directory
  if(chdir(directory_) != 0){
    carp(CARP_WARNING, "Cannot enter index directory %s", directory_);
    return false;
  }

  // create README file, with parameter informations
  FILE* readme = fopen("README", "w");
  writeReadmeFile(readme);
  fclose(readme);

  // create the index map & info
  FILE* info_out = fopen("crux_index_map", "w");
  writeHeader(info_out);
               
  // create text file of peptides.
  FILE* text_file = NULL;
  if (create_text_file) {
    carp(CARP_DEBUG, "Creating peptides.txt.");
    text_file = fopen("peptides.txt", "w");
  }
  
  // first, index target database
  index_database(database_, Index::index_file_prefix, info_out, text_file);
  fclose(info_out);
  info_out = NULL;
  if( text_file){
    fclose(text_file);
    text_file = fopen("decoy-peptides.txt", "w");
  }
  
  // now index decoy database without the info files
  info_out = fopen("crux_decoy_index_map", "w");
  writeHeader(info_out);
  index_database(decoy_database_, Index::decoy_index_file_prefix, 
                 info_out, text_file);//NULL);
  fclose(info_out);
  info_out = NULL;
  if( text_file){
    fclose(text_file);
    text_file = NULL;
  }
  
  //move out of index dir
  if( chdir("..") == -1 ){ 
    return false;
  }

  on_disk_ = true;
  return true;
}

/**
 * The steps of creating an index that are repeated for the target and
 * decoy databases.  For the decoy, the info_file and text_file can be NULL.
 */
void Index::index_database(
  Database* database, ///< the database to index
  const char* file_prefix, ///< name for the index files
  FILE* info_out, ///< index map
  FILE* text_file) ///< optional peptides file
{
  if( database == NULL ){
    return;
  }

  // get number of bins/files needed
  int* mass_limits = (int*)mycalloc(2, sizeof(int));
  long num_bins = getNumBinsNeeded(mass_limits);

#ifdef _MSC_VER
  if (num_bins > WIN_MAX_OPEN_FILES) {
    carp(CARP_FATAL, "Needed to open %ld index files. Windows allows a maximum of %ld open files\n", num_bins, WIN_MAX_OPEN_FILES); 
  }
#endif
  
  // create file handle array
  FILE** file_array = (FILE**)mycalloc(num_bins, sizeof(FILE*));

  // peptide array to store the peptides before serializing them all together
  Peptide*** peptide_array = (Peptide***)
    mycalloc(num_bins, sizeof(Peptide**));

  for(int sub_indx = 0; sub_indx < num_bins; ++sub_indx){
    peptide_array[sub_indx] = (Peptide**)
      mycalloc(MAX_PROTEIN_IN_BIN, sizeof(Peptide*));
  }
  // int array that stores the peptide count for each peptide array branch
  // this is used to determine when to output all peptides in buffer
  // does not represent the total peptide count in the bin
  // total count of peptide in bin is sotred in peptide_count_array
  int* bin_count = (int*)mycalloc(num_bins, sizeof(int));

  // create array that stores total count of peptides in each bin
  unsigned int* peptide_count_array = 
    (unsigned int*)mycalloc(num_bins, sizeof(unsigned int));

  // create database peptide_iterator
  DatabasePeptideIterator* peptide_iterator =
    new DatabasePeptideIterator(database, disk_constraint_, 
                                false);// don't parse all pep into memory

  long int file_idx = 0;
  int low_mass = mass_limits[0];
  long int count_peptide = 0;
  int mod_me = 1000;
  
  // iterate through all peptides
  while(peptide_iterator->hasNext()){    
    ++count_peptide;
    if(count_peptide % mod_me == 0){
      if( (count_peptide/10 ) == mod_me ){
        mod_me = mod_me * 10;
      }
      carp(CARP_INFO, "Reached peptide %d", (int)count_peptide);
    }

    Peptide* working_peptide = peptide_iterator->next();
    FLOAT_T working_mass = working_peptide->getPeptideMass();
    file_idx = (long int)((working_mass - low_mass) / mass_range_);

    // check if first time using this bin, if so create new file handle
    if(file_array[file_idx] == NULL){
      if(!generate_one_file_handle(file_array, file_idx, file_prefix)){
        fcloseall();
        carp(CARP_FATAL, 
             "Could not open file handle for index %d.", file_idx);
      }
    }
    
    // increment total peptide count of bin
    ++peptide_count_array[file_idx];

    // dump peptide in bin or temporary matrix
    dump_peptide(file_array, file_idx, working_peptide, 
                 peptide_array[file_idx], bin_count); 
  }

  carp(CARP_INFO, "Printing index");
  // dump all the left over peptides
  dump_peptide_all(file_array, peptide_array, bin_count, num_bins);

  // sort each bin  
  carp(CARP_INFO, "Sorting index");
  long bin_idx;
  for(bin_idx = 0; bin_idx < num_bins; ++bin_idx){
    carp(CARP_DETAILED_DEBUG, "Sorting bin %d", bin_idx);
    if(file_array[bin_idx] == NULL){
      continue;
    }
    // sort bin
    if((file_array[bin_idx] = sortBin(file_array[bin_idx], 
                                      bin_idx,  
                                      peptide_count_array[bin_idx], 
                                      text_file,
                                      file_prefix,
                                      database)) == NULL){
      fcloseall();
      carp(CARP_FATAL, "Failed to sort bin %i", bin_idx);
    }

    // print to crux_map
    if( info_out ){
      char* filename = get_crux_filename(bin_idx, file_prefix); 
      fprintf(info_out, "%s\t%.2f\t", filename, 
              mass_limits[0] + (bin_idx * mass_range_));
      fprintf(info_out, "%.2f\n", mass_range_);
      
      // free up heap
      std::free(filename);
    }
    fclose(file_array[bin_idx]);
  }
  
  std::free(mass_limits);
  std::free(file_array);
  std::free(peptide_count_array);
  delete (peptide_iterator);
}

/*********************************************
 * set and get methods for the object fields
 *********************************************/

/**
 *\returns A const pointer to the directory of the index
 */
const char* Index::getDirectory()
{
  return directory_;
}

/**
 * sets the directory of the index
 * index->directory must been initiailized
 */
void Index::setDirectory(
  const char* directory ///< the directory to add -in
  )
{
  std::free(directory_);
  directory_ = my_copy_string(directory);
}

/**
 *\returns a pointer to the database
 */
Database* Index::getDatabase(
  bool is_decoy)///< if true, return decoy database
{
  if( is_decoy ){
    return decoy_database_;
  }
  return database_;
}

/**
 *\returns a pointer to the peptides constraint
 */
PeptideConstraint* Index::getSearchConstraint() {

  return search_constraint_;
}

/**
 * \returns the disk constraint pointer from index
 */
PeptideConstraint* Index::getDiskConstraint() {
  return disk_constraint_;
}

/**
 * \brief Sets the peptide search constraint to be used by the
 * generate_peptides_iterator.  Makes a copy of the constraint pointer.
 * Deletes any existing search constraint. 
 */
void Index::setSearchConstraint(
  PeptideConstraint* constraint ///< Constraint for the next iterator
  )
{
  if( search_constraint_ ){
    PeptideConstraint::free(search_constraint_);
  }
  search_constraint_ = PeptideConstraint::copyPtr(constraint);
  // check that the new mass window is within the index
  double search_min = constraint->getMinMass();
  double search_max = constraint->getMaxMass();
  double index_min = disk_constraint_->getMinMass();
  double index_max = disk_constraint_->getMaxMass();
  
  if( search_min < index_min ){
    carp_once(CARP_WARNING, 
              "Minimum mass in the search range (%g) is below the index minimum (%g).",
              search_min, index_min);
    carp_once(CARP_WARNING, "This warning will not be repeated.")
  }
  if( search_max > index_max ){
    carp_once(CARP_WARNING, 
              "Maximum mass in the search range (%g) is above the index maximum (%g).",
              search_max, index_max);
    carp_once(CARP_WARNING, "This warning will not be repeated.")
  }
}

/**
 *\returns true if only allow unique peptides else false
 */
bool Index::getIsUnique()
{
  return is_unique_;
}

/**
 * \returns The type of decoys stored in the database.
 */
DECOY_TYPE_T Index::getDecoyType(){
  return decoys_;
}

/**
 * \returns the IndexMap for the target or decoy index
 */
IndexMap* Index::getIndexMap(
  bool decoy ///< return decoy index?
  ) {

  if (decoy) {
    if (decoy_index_map_ == NULL) {
      decoy_index_map_ = new IndexMap(this, true);
    }
    return decoy_index_map_;
  } else {
    if (index_map_ == NULL) {
      index_map_ = new IndexMap(this, false);
    }
    return index_map_;
  }


}


/****************************
 * bin_peptide_iterator
 ******************************/

/**
 * reads in the first peptide from the bin
 * parse the peptide, then adds it to the bin-peptide-iterator to return
 * \returns true if successfully initializes the bin_peptide_iterator
 */
bool initialize_bin_peptide_iterator(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< working bin_peptide_iterator -in/out
  )
{
  // BUG: This no longer distinguishes between eof and error in
  // parsing.  One fix would be for parse_peptide to return error code

  FILE* file = bin_peptide_iterator->index_file;
  Database* database = bin_peptide_iterator->database;// not decoy

  // allocate peptide to used to parse
  Peptide* peptide = Peptide::parse(file, database);

  if( peptide == NULL ){
    bin_peptide_iterator->peptide = NULL;
    bin_peptide_iterator->has_next = false;
    return true;
  }
  // set file pointer
  bin_peptide_iterator->index_file = file;
  bin_peptide_iterator->peptide = peptide;
  bin_peptide_iterator->has_next = true;
  return true;
}


/**
 * Instantiates a new bin_peptide_iterator from a gvien bin file handle.
 * \returns a new heap allocated bin_peptide_iterator object
 */
BIN_PEPTIDE_ITERATOR_T* new_bin_peptide_iterator(
  Index* index, ///< The index object which we are iterating over -in
  FILE* file, ///< the bin to parse peptides
  Database* database ///< the index's database to use (target/decoy)
  )
{
  // allocate a new index_peptide_iterator object
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator =
    (BIN_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(BIN_PEPTIDE_ITERATOR_T));
  
  // set index, file
  bin_peptide_iterator->index = index;
  bin_peptide_iterator->database = database;
  bin_peptide_iterator->index_file = file;
    
  if(!initialize_bin_peptide_iterator(bin_peptide_iterator)){
    carp(CARP_WARNING, "failed to initalize bin peptide iterator");
    bin_peptide_iterator->has_next = false;
  }
    
  return bin_peptide_iterator;
}


/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
Peptide* bin_peptide_iterator_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the bin_peptide_iterator to get peptide -in
  )
{
  Peptide* peptide_to_return = bin_peptide_iterator->peptide;

  // check if there's actually a peptide to return
  if(!bin_peptide_iterator_has_next(bin_peptide_iterator) ||
     bin_peptide_iterator->peptide == NULL){
    carp(CARP_FATAL, "bin_peptide_iterator, no peptides to return");
  }
  
  // setup the interator for the next peptide, if avaliable
  if(!initialize_bin_peptide_iterator(bin_peptide_iterator)){
    //die("failed to setup bin_peptide_iterator for next iteration");
    carp(CARP_WARNING, "am I fucking things up with the iterator???");
  }
 
  return peptide_to_return;
}

/**
 * The basic iterator functions.
 * check to see if the bin_peptide_iterator has more peptides to return
 *\returns true if there are additional peptides to iterate over, false if not.
 */
bool bin_peptide_iterator_has_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the bin_peptide_iterator to initialize -in
  )
{
  return bin_peptide_iterator->has_next; 
}

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_peptide_iterator(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the iterator to free -in
  )
{
  
  // if did not iterate over all peptides, free the last peptide not returned
  if(bin_peptide_iterator_has_next(bin_peptide_iterator)){
    delete bin_peptide_iterator->peptide;
  }

  free(bin_peptide_iterator);
}

/******************************
 * bin_sorted_peptide_iterator
 ******************************/

/**
 * Instantiates a new sorted_bin_peptide_iterator from a gvien bin
 * file handle. 
 * \returns a new heap allocated sorted_bin_peptide_iterator object
 */
BIN_SORTED_PEPTIDE_ITERATOR_T* new_bin_sorted_peptide_iterator(
  Index* index, ///< The index object which we are iterating over -in
  FILE* file, ///< the working bin file handler -in
  unsigned int peptide_count, ///< the total peptide count in the bin -in
  Database* database ///< the index's database to use (target/decoy) -in
  )
{
  // set peptide implementation to array peptide_src
  // this determines which peptide free method to use
  //  Peptide::setPeptideSrcImplementation(false);
  
  // create database sorted peptide iterator
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator =
    (BIN_SORTED_PEPTIDE_ITERATOR_T*)
    mycalloc(1, sizeof(BIN_SORTED_PEPTIDE_ITERATOR_T));

  // reset file to start
  rewind(file);

  // create bin_peptide_iterator
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator =
    new_bin_peptide_iterator(index, file, database);

  // create a sorted peptide iterator that will sort all the peptides 
  // from bin peptide_iterator
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator = 
    new_sorted_peptide_iterator_bin(bin_peptide_iterator, SORT_MASS, 
        index->getIsUnique(), peptide_count);

  // set sorted_peptide_iterator
  bin_sorted_peptide_iterator->sorted_peptide_iterator = 
    sorted_peptide_iterator;
  
  free_bin_peptide_iterator(bin_peptide_iterator); 
  // CHECK ME might wanna check this...

  return bin_sorted_peptide_iterator;
}

/**
 * The basic iterator functions.
 * \returns The next peptide in the index.
 */
Peptide* bin_sorted_peptide_iterator_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator 
    ///< the bin_peptide_iterator to get peptide -in
  )
{
  Peptide* peptide = sorted_peptide_iterator_next(bin_sorted_peptide_iterator->sorted_peptide_iterator);
  return peptide;
}

/**
 * The basic iterator functions.
 * check to see if the bin_sorted_peptide_iterator has more peptides to return
 *\returns true if there are additional peptides to iterate over, false if not.
 */
bool bin_sorted_peptide_iterator_has_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator ///< the bin_peptide_iterator to initialize -in
  )
{
  return sorted_peptide_iterator_has_next(bin_sorted_peptide_iterator->sorted_peptide_iterator);
}

/**
 * Frees an allocated bin_peptide_iterator object.
 */
void free_bin_sorted_peptide_iterator(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator ///< the iterator to free -in
  )
{
  free_sorted_peptide_iterator(bin_sorted_peptide_iterator->sorted_peptide_iterator);
  free(bin_sorted_peptide_iterator);
}




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
