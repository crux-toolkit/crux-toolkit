/************************************************************************//**
 * \file index.cpp
 * \brief Object for representing an index of a database
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>

#include "utils.h"
#include "crux-utils.h"

#include "peptide.h"
#include "Protein.h"
#include "index.h"
#include "carp.h"
#include "sorter.h"
#include "objects.h"
#include "peptide_constraint.h"
#include "database.h"
#include "protein_index.h"
#include "parameter.h"

// maximum proteins the index can handle
static const int MAX_PROTEIN = 30000;
static const int MAX_FILE_NAME_LENGTH = 300;
static const int NUM_CHECK_LINES = 8;
static const int MAX_PROTEIN_IN_BIN = 2500;
static const int MAX_FILE_SIZE_TO_USE_LIGHT_PROTEIN = 500000000;
static const int MAX_PARSE_COUNT = 3;
static const int SLEEP_DURATION = 5;

#ifdef DARWIN
// OS X doesn't support fcloseall()
// Consider tracking open file handles
// rather then simply relying on process
// termination to close files.
#define fcloseall()
#endif

// global variable to store the temp directory
// used for deleting directory when SIGINT
char temp_folder_name[12] = "";

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

/* Private Functions */
// BOOLEAN_T set_index_fields_from_disk(INDEX_T* index);
void set_index_field_from_map(INDEX_T* index, char* line);
BOOLEAN_T check_index_constraints(INDEX_T* index);

/* Private Data Types */
/**
 * \struct index
 * \brief A index of a database
 *
 * An index consists of three parts: a file with protein sequences, 
 * a set of files that point into the sequence file, and a map file
 * listing which files contain pointers for which peptide masses.  All
 * files are all contained in a single directory.
 *
 * The directory is also the name of the index given by the user.
 * There are no assumptions made about the name of the index
 * (directory) and the name of the fasta file from which it was
 * created.  The binary version of the fasta is named [fasta
 * filename]-binary-fasta.  The file with a list of index files and
 * the masses they contain is crux_index_map.
 *
 * The index has two peptide constraints: one to define ALL the
 * peptides that are indexed and one to define the peptides being
 * searched for.  When an index is being created, the search
 * constraint is NULL.  When an index is being loaded from disk, the
 * disk_constraint is set according to the values in the header of
 * crux_index_map.  Each time a peptide generator is created for the
 * index, the search_constraint is set (any existing ones being
 * deleted).  
 *
 * Bounds checking for a search is done at two points.  When the index
 * is loaded from disk, the global constraints for the search are
 * compared to the constraints for the index and the process is
 * terminated if they are not compatible. (see
 * check_index_constraints)  When a peptide iterator is created and a 
 * new search_constraint is set, only the min and max mass is checked
 * (as no other parts of the constraint change in the course of a
 * search).  If the iterator requests a mass range partly outside the
 * bounds of the index, a warning is printed (one for over max and one
 * for under min).
 */
struct index{
  int num_pointers; ///< The number of pointers to this index.
  DATABASE_T* database; ///< The database that has been indexed.
  char* directory; ///< The directory containing the indexed files
  PEPTIDE_CONSTRAINT_T* disk_constraint;///< Defines peptides on disk
  PEPTIDE_CONSTRAINT_T* search_constraint;///< Defines peptides being searched
  BOOLEAN_T on_disk; ///< Does this index exist on disk yet?
  FLOAT_T mass_range;  ///< the range of masses contained in each index file -in
  BOOLEAN_T is_unique; ///< only unique peptides? -in
};    

/**
 * \struct index_file
 * \brief A struct that contains the information of each file
 */
struct index_file{
  char* filename;  ///< The file name that contain the peptides
  FLOAT_T start_mass; ///< the start mass limit in this file
  FLOAT_T interval;   ///< the interval of the peptides in this file
};
typedef struct index_file INDEX_FILE_T;


/**
 * \struct index_peptide_iterator
 * \brief An iterator to iterate over the peptides in a database
 */
struct index_peptide_iterator{
  INDEX_T* index; ///< The index object which we are iterating over
  INDEX_FILE_T* index_files[MAX_INDEX_FILES]; 
  ///< the index file array that contain information of each index file 
  int total_index_files; ///< the total count of index_files
  int current_index_file; ///< the index file open or one to open next 
  FILE* index_file; ///< The current file stream that we are reading from
  BOOLEAN_T has_next; ///< Is there another peptide?
  PEPTIDE_T* peptide; ///< the next peptide to return
};    

/**
 * \struct index_filtered_peptide_iterator
 * \brief An iterator to filter out the peptides wanted from the
 * index_peptide_iterator 
 */
struct index_filtered_peptide_iterator{
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator;///< Core peptide iterator
  BOOLEAN_T has_next; ///< Is there another peptide?
  PEPTIDE_T* peptide; ///< the next peptide to return
};    

/**
 * \struct bin_peptide_iterator
 * \brief An iterator to iterate over the peptides in a bin (one file
 * handle) 
 */
struct bin_peptide_iterator{
  INDEX_T* index; ///< The index object which we are iterating over
  FILE* index_file; ///< The current file stream that we are reading from
  BOOLEAN_T has_next; ///< Is there another peptide?
  PEPTIDE_T* peptide; ///< the next peptide to return
  BOOLEAN_T use_array; 
  ///< Use array peptide_src or link list peptide_src when parsing peptides
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
 * Public functions for INDEX_T
 ************/

//TODO: shouldn't allocate be private?
/**
 * \returns An (empty) index object.
 */
INDEX_T* allocate_index(void){
  INDEX_T* index = (INDEX_T*)mycalloc(1, sizeof(INDEX_T));
  index->num_pointers = 1;
  return index;
}

/**
 * Helper function for scandir to find files with names ending in
 * "-binary-fasta"
 * \returns 1 if filename ends in -binary-fasta, else 0 
 */
#ifdef DARWIN
int is_binary_fasta_name(struct dirent *entry){
#else
int is_binary_fasta_name(const struct dirent *entry){
#endif


  const char* filename = entry->d_name; //w/o const gcc warning
  const char* suffix = "-binary-fasta";

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

// TODO (BF 27-Feb-08): find what I used to generate the name and put
// the two functions near each other.
/**
 * \brief Looks in given directory for a file ending in
 * "-binary-fasta" and returns a heap-allocated string of the full
 * name including the index directory.
 *
 * Exits with error if index->directory does not exist, no file
 * *-binary-fasta exists, or more than one *binary-fasta file exists.
 * \returns A string with the name of the existing binary fasta file
 * for this index. 
 */
char* get_index_binary_fasta_name(const char* index_name){
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

  free(namelist[0]);
  free(namelist);
  free(fasta_name);

  return fasta_name_path;

}

/**
 * \brief Initializes a new INDEX_T object by setting its member
 * variables to the values in the index on disk.
 *
 * Assumes that the name of the index (directory) has been set by
 * caller.  Checks that the directory exists.  Reads the mapfile
 * therein and gets the constraint values from the header.  Returns
 * TRUE if all of the above are true, else prints an error and returns
 * FALSE. 
 *
 * \returns TRUE if fields were successfully set, FALSE if index on
 * disk could not be used to initialize the object.
 */
BOOLEAN_T set_index_fields_from_disk(
  INDEX_T* index  ///< Index to set -out                       
  )
{
  assert(index->directory != NULL );
  DIR* check_dir = opendir(index->directory);
  
  if(check_dir == NULL){
    carp(CARP_ERROR, "Failed to access index directory %s", index->directory);
    return FALSE;
  }else{
    closedir(check_dir);
  }

  // create constraint
  index->disk_constraint = allocate_peptide_constraint();
  index->search_constraint = NULL;

  // open map file
  char* map_filename = get_full_filename(index->directory, "crux_index_map");
  carp(CARP_DEBUG, "Opening index map file '%s'", map_filename);
  FILE* map_file = fopen(map_filename, "r");
  if(map_file == NULL){
    carp(CARP_ERROR, "Could not open index map file %s", map_filename);
    return FALSE;
  }

  // read header lines (beginning with #)
  char* line = NULL;
  size_t buf_length = 0;
  int line_length = getline(&line, &buf_length, map_file);
  while( line_length != -1 && line[0] == '#' ){
    set_index_field_from_map(index, line);
    line_length = getline(&line, &buf_length, map_file);
  }
  fclose(map_file);
  free(map_filename);
  myfree(line);

  index->on_disk = TRUE;
  return TRUE;
}

/**
 * \brief Private function to take a header line from index_map and
 * set the appropriate value in the index.
 */
void set_index_field_from_map(INDEX_T* index, char* line){
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
    set_peptide_constraint_min_mass(index->disk_constraint, value);
  }
  else if(strcmp("max_mass:", trait_name) == 0){
    set_peptide_constraint_max_mass(index->disk_constraint, value);
  }
  else if(strcmp("min_length:", trait_name) == 0){
    set_peptide_constraint_min_length(index->disk_constraint, (int)value);
  }
  else if(strcmp("max_length:", trait_name) == 0){
    set_peptide_constraint_max_length(index->disk_constraint, (int)value);
  }
  else if(strcmp("peptide_type:", trait_name) == 0){
    //    set_peptide_constraint_peptide_type(index->disk_constraint, value);
    carp(CARP_FATAL, "This index uses the obsolete 'peptide_type'. "
         "Rebuild with current version of crux to use 'enzyme_type'.");
  }
  else if(strcmp("enzyme_type:", trait_name) == 0){
    set_peptide_constraint_enzyme(index->disk_constraint, (ENZYME_T)value);
  }
  else if(strcmp("digest_type:", trait_name) == 0){
    set_peptide_constraint_digest(index->disk_constraint, (DIGEST_T)value);
  }
  else if(strcmp("missed_cleavages:", trait_name) == 0){
   set_peptide_constraint_num_mis_cleavage(index->disk_constraint, (int)value);
  }
  else if(strcmp("mass_type:", trait_name) == 0){
   set_peptide_constraint_mass_type(index->disk_constraint, (MASS_TYPE_T)value);
  }
  else if(strcmp("unique_peptides:", trait_name) == 0){
    index->is_unique = (BOOLEAN_T)value;
  }
  else if(strcmp("target_mass_range_for_index_file:", trait_name) == 0){
    index->mass_range = value;
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
 * \returns TRUE if all constraints will work with those in parameter.c.
 */
BOOLEAN_T check_index_constraints(INDEX_T* index){
  double min_mass = get_peptide_constraint_min_mass(index->disk_constraint);
  double max_mass = get_peptide_constraint_max_mass(index->disk_constraint);
  double min_len = get_peptide_constraint_min_length(index->disk_constraint);
  double max_len = get_peptide_constraint_max_length(index->disk_constraint);
  int missed_cleavages = 
    get_peptide_constraint_num_mis_cleavage(index->disk_constraint);
  //  PEPTIDE_TYPE_T cleavage_type = 
  //    get_peptide_constraint_peptide_type(index->disk_constraint);
  ENZYME_T enzyme = get_peptide_constraint_enzyme(index->disk_constraint);
  DIGEST_T digestion = get_peptide_constraint_digest(index->disk_constraint);
  MASS_TYPE_T mass_type = get_peptide_constraint_mass_type(index->disk_constraint);
  //  BOOLEAN_T unique = index->is_unique;

  BOOLEAN_T success = TRUE;
  const char* param;
  if(min_mass > get_double_parameter("min-mass")){
    success = FALSE;
    param = "min-mass";
  }else if(max_mass < get_double_parameter("max-mass")){
    success = FALSE;
    param = "max-mass";
  }else if(min_len > get_int_parameter("min-length")){
    success = FALSE;
    param = "min-length";
  }else if(max_len < get_int_parameter("max-length")){
    success = FALSE;
    param = "max-length";
  }else if(missed_cleavages < get_int_parameter("missed-cleavages")){
    success = FALSE;
    param = "missed-cleavages";
  }else if(mass_type != get_mass_type_parameter("isotopic-mass")){
    success = FALSE;
    param = "isotopic-mass";
    /*
  }else if(unique != get_boolean_parameter("unique-peptides")){
    // TODO (BF 07-24-08): would like to change so that non-unique
    // index could return unique peptides in generate-peptides
    // only if index is unique(1) and requesting not-unique (0) is a problem
    success = FALSE;
    param = "unique-peptides";
    */
    /*
  }else if(cleavage_type < get_peptide_type_parameter("cleavages")){
    // tryptic < partial < any  BUG for N- or C-tryptic)
    success = FALSE;
    param = "cleavages";
  }
    */
  }else if(enzyme != get_enzyme_type_parameter("enzyme")){
    printf("constraint e: %i, param e: %i.", (int)enzyme, (int)get_enzyme_type_parameter("enzyme"));
    success = FALSE;
    param = "enzyme";
  }else if(digestion < get_digest_type_parameter("digestion")){
    // full < partial < non-specific
    success = FALSE;
    param = "digestion";
  }

  if(success == FALSE){
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
void set_index_fields(
  INDEX_T* index,  ///< Index to set -out                       
  const char* output_dir,      ///< The name of the new index
  PEPTIDE_CONSTRAINT_T* constraint,  
  ///< Constraint which these peptides satisfy -in
  FLOAT_T mass_range,  
  ///< the range of mass that each index file should be partitioned into -in
  BOOLEAN_T is_unique ///< only unique peptides? -in
  )
{
  carp(CARP_DEBUG, "Setting index fields");
  DIR* check_dir = NULL;
  
  if((check_dir = opendir(output_dir)) != NULL){
    set_index_on_disk(index, TRUE);
    closedir(check_dir);
  }
  else{
    set_index_on_disk(index, FALSE);
  }

  carp(CARP_DEBUG, "Index on disk is '%i' for dir '%s'", 
       (int)index->on_disk, output_dir);

  // set each field
  set_index_directory(index, output_dir);
  index->disk_constraint = constraint;
  index->search_constraint = NULL;
  //index->on_disk = FALSE; // this breaks overwrite of create index
  index->mass_range = mass_range;  
  index->is_unique = is_unique;

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
INDEX_T* new_index(
  const char* fasta_filename, ///< The fasta file
  const char* output_dir,     ///< The name of the new index
  PEPTIDE_CONSTRAINT_T* constraint,  
    ///< Constraint which these peptides will satisfy
  FLOAT_T mass_range  
    ///< the range of mass that each index file should be partitioned into
  )
{
  carp(CARP_DETAILED_DEBUG, "Creating new index to be named %s", output_dir);

  INDEX_T* index = allocate_index();
  DATABASE_T* database = NULL;
  
  // Initially, create a database that does not use memory mapping (FALSE)
  // Once binary fasta file has been creaated this will change to a
  // memory mapped database
  database = new_database(fasta_filename, FALSE);

  // set database, has not been parsed
  set_index_database(index, database);
  BOOLEAN_T is_unique = TRUE;
  set_index_fields(index, output_dir,
                   constraint, mass_range, is_unique);

  return index;
}         

/**
 * \brief Create an index object to represent an index already on disk.
 * Used by generate-peptides and search-for-matches.  If the proper
 * files are not found in the directory, warns and returns NULL.
 *
 * \returns A new index object ready for search.
 */
INDEX_T* new_index_from_disk(
  const char* index_name  ///< The directory containing the index
  )
{
  // find the database file, open it and point to it
  // check that it works with global constraints

  // allocate index and name it
  INDEX_T* search_index = allocate_index();
  set_index_directory(search_index, index_name);
  
  // this should find the map file and get constraint values
  if(!set_index_fields_from_disk(search_index)){
    carp(CARP_FATAL, "Could not create index %s from disk", index_name);
  }

  // check that the index constraints are OK for this run
  if(!check_index_constraints(search_index)){
    carp(CARP_FATAL, "Index cannot produce the requested peptides. " \
     "Change the search parameters to match the index or create new index.");
  }

  // get binary fasta file name with path to crux directory 
  char* binary_fasta = get_index_binary_fasta_name(index_name);
  
  // check if input file exist
  if(access(binary_fasta, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist for crux index.", 
        binary_fasta);
  }
  
  // now create a database, using binary fasta file
  search_index->database = new_database(binary_fasta, TRUE);
  
  if(!parse_database(search_index->database)){
    carp(CARP_FATAL, "Failed to parse database, cannot create new index");
  }

  // free string
  free(binary_fasta);
  
  return search_index;
}

int get_index_num_proteins(INDEX_T* index){

  return get_database_num_proteins(index->database);
}

/**
 * Frees an allocated index object.
 */
int get_index_pointer_count(
  INDEX_T* index
  ){
  return index->num_pointers;
}


/**
 * Frees an allocated index object.
 */
INDEX_T* copy_index_ptr(
  INDEX_T* index
  ){
  index->num_pointers++;
  return index;
}

/**
 * Frees an allocated index object, but only if pointer count is equal to 1.
 */
void free_index(
  INDEX_T* index
  )
{
  if(index == NULL ){
    return;
  }

  if (index->num_pointers > 1){
    index->num_pointers--;
  } else {
    carp(CARP_DEBUG, "Freeing index");
    if (index->database != NULL){
      carp(CARP_DEBUG, "Freeing index database");
      free_database(index->database);
    }
    if (index->disk_constraint != NULL){
      carp(CARP_DEBUG, "Freeing index disk constraint");
      free_peptide_constraint(index->disk_constraint);
    }
    if (index->search_constraint != NULL){
      carp(CARP_DEBUG, "Freeing index search constraint");
      free_peptide_constraint(index->search_constraint);
    }
    free(index->directory);
    free(index);
  }
}

/**
 * write to the file stream various information of the
 * index files created
 */
BOOLEAN_T write_header(
  INDEX_T* index, ///< the working index -in
  FILE* file ///< out put stream for crux_index_map -in
  )
{
  time_t hold_time;
  hold_time = time(0);
  PEPTIDE_CONSTRAINT_T* constraint = index->disk_constraint;
  
  
  fprintf(file, "#\tmin_mass: %.2f\n", get_peptide_constraint_min_mass(constraint));
  fprintf(file, "#\tmax_mass: %.2f\n", get_peptide_constraint_max_mass(constraint));
  fprintf(file, "#\tmin_length: %d\n", get_peptide_constraint_min_length(constraint));
  fprintf(file, "#\tmax_length: %d\n", get_peptide_constraint_max_length(constraint));
  //fprintf(file, "#\tpeptide_type: %d\n", get_peptide_constraint_peptide_type(constraint));
  fprintf(file, "#\tenzyme_type: %d\n", get_peptide_constraint_enzyme(constraint));
  fprintf(file, "#\tdigest_type: %d\n", get_peptide_constraint_digest(constraint));
  fprintf(file, "#\tmissed_cleavages: %d\n", get_peptide_constraint_num_mis_cleavage(constraint));
  fprintf(file, "#\tmass_type: %d\n", get_peptide_constraint_mass_type(constraint));
  fprintf(file, "#\tunique_peptides: %d\n", get_index_is_unique(index));
  
  fprintf(file, "#\tCRUX index directory: %s\n", index->directory);
  fprintf(file, "#\ttime created: %s",  ctime(&hold_time)); 
  fprintf(file, "#\ttarget_mass_range_for_index_file: %.2f\n", index->mass_range);
  
  return TRUE;
}

/**
 * write to the file stream various information of the
 * index files created in human readable format
 *\returns TRUE if successfully creates README file, else FALSE
 */
BOOLEAN_T write_readme_file(
  INDEX_T* index, ///< the working index -in
  FILE* file ///< out put stream for README file -in
  )
{
  time_t hold_time;
  hold_time = time(0);
  PEPTIDE_CONSTRAINT_T* constraint = index->disk_constraint;
  char* fasta_file = get_database_filename(index->database);
  char* fasta_file_no_path = parse_filename(fasta_file);
  
  fprintf(file, "#\ttime created: %s",  ctime(&hold_time)); 
  fprintf(file, "#\tfasta file: %s\n",  fasta_file_no_path); 
  fprintf(file, "#\tmin_mass: %.2f\n",
          get_peptide_constraint_min_mass(constraint));
  fprintf(file, "#\tmax_mass: %.2f\n",
          get_peptide_constraint_max_mass(constraint));
  fprintf(file, "#\tmin_length: %d\n",
          get_peptide_constraint_min_length(constraint));
  fprintf(file, "#\tmax_length: %d\n",
          get_peptide_constraint_max_length(constraint));

  /*
  PEPTIDE_TYPE_T type = get_peptide_constraint_peptide_type(constraint);
  char type_str[64];
  peptide_type_to_string(type, type_str);
  fprintf(file, "#\tpeptide_type: %s\n", type_str);
  */
  char* value_str = 
    enzyme_type_to_string(get_peptide_constraint_enzyme(constraint));
  fprintf(file, "#\tenzyme_type: %s\n", value_str);
  free(value_str);
  value_str = 
    digest_type_to_string(get_peptide_constraint_digest(constraint));
  free(value_str);

  fprintf(file, "#\tmissed_cleavage: %s\n", 
        (get_peptide_constraint_num_mis_cleavage(constraint)? "true":"false"));
  fprintf(file, "#\tmass_type: %s\n", 
          (get_peptide_constraint_mass_type(constraint)==AVERAGE? 
           "average":"mono"));
  fprintf(file, "#\tunique peptides: %s\n", 
          (get_index_is_unique(index)? "unique":"redundant"));
  fprintf(file, "#\tCRUX index directory: %s\n", index->directory);
  fprintf(file, "#\ttarget mass range for index file: %.2f\n", 
          index->mass_range);
  
  free(fasta_file);
  free(fasta_file_no_path);

  return TRUE;
}


/**
 * heap allocated, users must free
 * \returns a temporary directory name template
 */
char* make_temp_dir_template(void){
  char* dir_template = (char*)mycalloc(12, sizeof(char));
  strcpy(dir_template, "crux_XXXXXX");
  return dir_template;
}

/**
 * heap allocated filename
 *\returns the filename for the given index
 */
char* get_crux_filename(
  long bin_idx,  ///< the bin_indx name you want -in
  int part  ///< the what sub part of the dir is it? only needed when spliting -in
  )
{
  char* file_num = 0;
  const char* filename_tag = "crux_index_";
  char* filename = NULL;

  // quiet compiler
  part = part;

  // add functionallity to make _1, _2, _3
  file_num = int_to_char(bin_idx + 1);
  filename = cat_string(filename_tag, file_num);

  free(file_num);
  return filename;
}

/**
 * Calculate the total number of bins( file handlers) that will be needed
 * \returns the total number of bins needed
 */
long get_num_bins_needed(
  INDEX_T* index, ///< working index -in
  int* mass_limits  ///< an array that holds the min/max mass limit -in
  )
{
  int min_length = get_peptide_constraint_min_length(index->disk_constraint);
  int max_length = get_peptide_constraint_max_length(index->disk_constraint);
  FLOAT_T min_mass = get_peptide_constraint_min_mass(index->disk_constraint);
  FLOAT_T max_mass = get_peptide_constraint_max_mass(index->disk_constraint);
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

  num_bins = (int)(((max_mass_limit - min_mass_limit) / index->mass_range) + 1);  // check..

  if (num_bins > MAX_INDEX_FILES){
    num_bins = MAX_INDEX_FILES;
  }
  return num_bins;
}                         

/**
 * user MUST set the unix system max allowed file handlers enough to allow this procsess
 * check and change on command line by "ulimit -n", need root permission to change...
 *generates all the file handlers(bins) that are needed
 *\returns TRUE, if successfully opened all needed bins, else FALSE
 */
BOOLEAN_T generate_file_handlers(
  FILE** file_array,  ///< the file handler array -out
  long num_bins  ///< total number of bins needed
  )
{
  long bin_indx = 0;
  FILE* file = NULL;
  char* filename = NULL;
  
  // create all the file handlers need for create index
  for(; bin_indx < num_bins; ++bin_indx){
    filename = get_crux_filename(bin_indx, 0);
    file = fopen(filename, "w+" );
    free(filename);
    if(file == NULL){
      carp(CARP_WARNING, "cannot open all file handlers needed");
      free(file_array);
      return FALSE;
    }
    (file_array)[bin_indx] = file;
  }
  
  return TRUE;
}

/**
 * user MUST set the unix system max allowed file handlers enough to allow this procsess
 * check and change on command line by "ulimit -n", need root permission to change...
 *generates all the file handlers(bins) that are needed
 *\returns TRUE, if successfully opened bin, else FALSE
 */
BOOLEAN_T generate_one_file_handler(
  FILE** file_array,  ///< the file handler array -out                            
  long bin_index ///< the bin index to create a file handler -in
  )
{
  FILE* file = NULL;
  char* filename = NULL;
  
  // create the file handler needed for create index
  filename = get_crux_filename(bin_index, 0);
  file = fopen(filename, "w+" );
  free(filename);
  if(file == NULL){
    carp(CARP_WARNING, "cannot open all file handlers needed");
    return FALSE;
  }
  file_array[bin_index] = file;
  return TRUE;
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
FILE* sort_bin(
  FILE* file, ///< the working file handle to the bin -in
  long bin_idx, ///< bin index in the file array -in
  INDEX_T* index, ///< working index -in
  unsigned int peptide_count, ///< the total peptide count in the bin -in
  FILE* text_file
  )
{
  char* filename = NULL;

  // check if file is empty
  if(peptide_count == 0){
    return file;
  }

  BIN_SORTED_PEPTIDE_ITERATOR_T* peptide_iterator =
    new_bin_sorted_peptide_iterator(index, file, peptide_count);
  PEPTIDE_T* working_peptide = NULL;
  
  // get the filename for this file bin
  filename = get_crux_filename(bin_idx, 0);
  // close unsorted bin
  fclose(file);
  // create new bin which will be sorted 
  file = fopen(filename, "w" );
  
  // serialize all peptides in sorted order
  while(bin_sorted_peptide_iterator_has_next(peptide_iterator)){
    working_peptide = bin_sorted_peptide_iterator_next(peptide_iterator);
    serialize_peptide(working_peptide, file, text_file);
    free_peptide(working_peptide);
  }
  
  free(filename);
  free_bin_sorted_peptide_iterator(peptide_iterator);

  return file;
}

/**
 * Stores the peptide in the correct bin.  If bin exceeds
 * MAX_PROTEIN_IN_BIN, then serialize all peptides in the bin.
 *
 * \returns TRUE if successful in storing the peptide or serializing
 * peptides, else FALSE .
 */
static BOOLEAN_T dump_peptide(
  FILE** file_array,   ///< the working file handler array to the bins -in/out
  long int file_idx,   ///< the index of the file the peptide belongs to -in
  PEPTIDE_T* working_peptide, ///< the peptide to be stored -in
  PEPTIDE_T** peptide_array, ///< hold peptides before they're serialized -out
  int* bin_count       ///< an array of the counts of peptides in each bin -in
  )
{  
  int peptide_idx = 0;
  FILE* file = NULL;
  int current_count;
  
  // if the peptide count is over the limit
  if((current_count = bin_count[file_idx]) > MAX_PROTEIN_IN_BIN){
    file = file_array[file_idx];
    // print out all peptides
    for(peptide_idx = 0; peptide_idx < current_count; ++peptide_idx){
      serialize_peptide(peptide_array[peptide_idx], file, NULL);
      free_peptide(peptide_array[peptide_idx]);
      // FIXME this should not be allowed
    }
    serialize_peptide(working_peptide, file, NULL);
    free_peptide(working_peptide);
    // FIXME this should not be allowed
    bin_count[file_idx] = 0;
  }
  // if the peptide count is bellow the limit
  else{
    // store peptide in peptide array , these peptides will be printed later togehter
    peptide_array[(bin_count[file_idx])] = working_peptide;
    ++bin_count[file_idx];
  }
  return TRUE;
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
 * \returns TRUE, if successful in serializing all peptides, else FALSE.
 */
static BOOLEAN_T dump_peptide_all(
  FILE** file_array,   ///< the working file handle array to the bins -out
  PEPTIDE_T*** peptide_array, ///< the array of pre-serialized peptides -in
  int* bin_count,      ///< the count array of peptides in each bin -in
  int num_bins         ///< the total number of bins -in
  )
{  
  int peptide_idx = 0;
  FILE* file = NULL;
  PEPTIDE_T** working_array = NULL;
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
      serialize_peptide(working_array[peptide_idx], file, NULL);
      free_peptide(working_array[peptide_idx]);
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
  return TRUE;
}

/***
 * This function does the following things...
 * 1. create binary fasta file in temporary directory
 * 2. transform database into memory mapped database from text base database
 * 3. then, parse database
 * Called while the cwd is the temp directory in which the index is
 * being made.
 *
 *\returns TRUE, if all processes are successful, else FALSE
 */
BOOLEAN_T transform_database_to_memmap_database(
  INDEX_T* index ///< An allocated index -in/out
  )
{
  char* binary_fasta = NULL; // BF not used after set

  // get the fasta file name with correct path
  char* fasta_file = cat_string("../", 
                              get_database_filename_pointer(index->database));

  // create binary fasta file inside temp directory
  if(!create_binary_fasta_in_cur(fasta_file,
                               get_database_filename_pointer(index->database),
                               &binary_fasta)){
    /*
    // remove directory
    chdir("..");
    clean_up(1);
    // free index
    free_index(index);
    free(fasta_file);
    carp(CARP_FATAL, "Failed to create protein index on disk");
    */
    carp(CARP_ERROR, "Failed to create protein index on disk");
    return FALSE;
  }
  
  // change name of file to binary fasta
  set_database_filename(index->database, fasta_file);

  // check if already parsed
  if(!get_database_is_parsed(index->database)){
    carp(CARP_DEBUG,
       "Database was not parsed after creating the binary fasta, parsing now");

    if(!parse_database(index->database)){
      /*
      carp(CARP_FATAL, "failed to parse database, cannot create new index");
      free(index);
      free(fasta_file);
      free(binary_fasta);
      fcloseall();
      exit(1);
      */
      carp(CARP_ERROR, "Failed to parse database, cannot create new index");
      return FALSE;
    }
  }
  
  // free file name
  free(fasta_file);
  free(binary_fasta);

  return TRUE;
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
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T create_index(
  INDEX_T* index, ///< An allocated index -in/out
  BOOLEAN_T create_text_file ///< Should an ASCII text file be create? -in
  )
{
  // the file stream where the index creation information is sent
  FILE* info_out = NULL; 
  FILE* readme = NULL;
  // new stuff
  char* temp_dir_name = NULL;
  FILE** file_array = NULL;
  int* mass_limits = (int*)mycalloc(2, sizeof(int));
  long num_bins = 0;
  DATABASE_PEPTIDE_ITERATOR_T* peptide_iterator = NULL;
  PEPTIDE_T* working_peptide = NULL;
  FLOAT_T working_mass;
  char* filename = NULL;
  FLOAT_T mass_range = index->mass_range;
  unsigned int* peptide_count_array = NULL;
  BOOLEAN_T replace_index = FALSE;

  carp(CARP_DEBUG, "Creating index");
  // check if already created index
  if(index->on_disk){
    if(get_boolean_parameter("overwrite")){
      replace_index = TRUE;
      carp(CARP_DEBUG, "Will be replacing existing index");
      // wait to delete until index is successfully created
    }else{ // this should have already been checked, but...
      carp(CARP_FATAL, "Index '%s' already exists.  " \
      "Use --overwrite T to replace", index->directory);
    }
    //change the ondisk status?
  }
  
  // create temporary directory
  // temp_dir_name = "foo"; // CYGWIN
  if(mkdir(temp_dir_name, S_IRWXO) != 0){
    if((temp_dir_name = mkdtemp(make_temp_dir_template()))== NULL){
      carp(CARP_WARNING, "Cannot create temporary directory");
      return FALSE;
    }
  }
  
  // copy temporary folder name for SIGINT cleanup purpose
  strncpy(temp_folder_name, temp_dir_name, 12); 

  if(! transform_database_text_to_memmap(index->database, 
                                         temp_dir_name) ){
    clean_up(1);
    carp(CARP_FATAL, "Failed to create binary database from text fasta file");
  }

  // move into temporary directory
  if(chdir(temp_dir_name) != 0){
    carp(CARP_WARNING, "Cannot enter temporary directory");
    return FALSE;
  }

  // 1. create binary fasta file in temporary directory
  // 2. transform database into memory mapped database from text base database
  // 3. then, parse database
  /*
  if(! transform_database_to_memmap_database(index)){
    clean_up(1);
    // free ??
    exit(1);
    } */
 
  // get number of bins needed
  num_bins = get_num_bins_needed(index, mass_limits);
  
  // create file handle array
  file_array = (FILE**)mycalloc(num_bins, sizeof(FILE*));

  // peptide array to store the peptides before serializing them all together
  PEPTIDE_T*** peptide_array = (PEPTIDE_T***)
    mycalloc(num_bins, sizeof(PEPTIDE_T**));
  int sub_indx;
  for(sub_indx = 0; sub_indx < num_bins; ++sub_indx){
    peptide_array[sub_indx] = (PEPTIDE_T**)
      mycalloc(MAX_PROTEIN_IN_BIN, sizeof(PEPTIDE_T*));
  }
  // int array that stores the peptide count for each peptide array branch
  // this is used to determine when to output all peptides in buffer
  // does not represent the total peptide count in the bin
  // total count of peptide in bin is sotred in peptide_count_array
  int* bin_count = (int*)mycalloc(num_bins, sizeof(int));

  // create array that stores total count of peptides in each bin
  peptide_count_array = 
    (unsigned int*)mycalloc(num_bins, sizeof(unsigned int));

  // create README file, with parameter informations
  readme = fopen("README", "w");
  write_readme_file(index, readme);
  fclose(readme);

  // create text file of peptides.
  FILE* text_file = NULL;
  if (create_text_file) {
    carp(CARP_DEBUG, "Creating peptides.txt.");
    text_file = fopen("peptides.txt", "w");
  }
  
  // create the index map & info
  info_out = fopen("crux_index_map", "w");
  write_header(index, info_out);
                    
  // create database peptide_iterator
  peptide_iterator =
    new_database_peptide_iterator(index->database, index->disk_constraint, 
                                  false);// don't parse all pep into memory

  long int file_idx = 0;
  int low_mass = mass_limits[0];
  long int count_peptide = 0;
  int mod_me = 1000;
  
  // iterate through all peptides
  while(database_peptide_iterator_has_next(peptide_iterator)){    
    ++count_peptide;
    if(count_peptide % mod_me == 0){
      if( (count_peptide/10 ) == mod_me ){
        mod_me = mod_me * 10;
      }
      carp(CARP_INFO, "Reached peptide %d", (int)count_peptide);
    }

    working_peptide = database_peptide_iterator_next(peptide_iterator);
    working_mass = get_peptide_peptide_mass(working_peptide);
    file_idx = (long int)((working_mass - low_mass) / mass_range);

    // check if first time using this bin, if so create new file handle
    if(file_array[file_idx] == NULL){
      if(!generate_one_file_handler(file_array, file_idx)){
        carp(CARP_ERROR, 
             "Exceeded filehandle limit on system with %d files", file_idx);
        fcloseall();
        return FALSE;
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
    if((file_array[bin_idx] = sort_bin(file_array[bin_idx], bin_idx, index, 
                                       peptide_count_array[bin_idx], 
                                       text_file)) == NULL){
      carp(CARP_WARNING, "Failed to sort bin %i", bin_idx);
      fcloseall();
      return FALSE;
    }

    // TODO add splitting of files if necessary
    // print to crux_map
    filename = get_crux_filename(bin_idx, 0); 
      ///< 0 can change if we need to split the file
    fprintf(info_out, "%s\t%.2f\t", filename, 
        mass_limits[0] + (bin_idx * index->mass_range));
    fprintf(info_out, "%.2f\n", index->mass_range);

    // free up heap
    free(filename);
    fclose(file_array[bin_idx]);
  }
  
  // close crux_index_map, free heap allocated objects
  if (create_text_file) {
    fclose(text_file);
  }
  fclose(info_out);
  free(mass_limits);
  free(file_array);
  free(peptide_count_array);
  free_database_peptide_iterator(peptide_iterator);

  if( chdir("..") == -1 ){ //move out of temp dir
    return FALSE;
  }

  // rename temporary direcotry to final directory name
  // if replacing an existing index, remove current files and delete
  // dir
  if( replace_index == TRUE ){
    carp(CARP_DEBUG,"About to delete existing index directory, %s", 
         index->directory);
    delete_dir(index->directory);
  }
  if(rename(temp_dir_name, index->directory) != 0){
    carp(CARP_WARNING, "Cannot rename directory");
    return FALSE;
  }

  free(temp_dir_name);

  // set permission for the directory
  chmod(index->directory, S_IRWXU+S_IRWXG+S_IROTH+S_IXOTH);

  index->on_disk = TRUE;
  return TRUE;
}

/**
 * Does this index exist on disk?
 *
 * \returns TRUE if it does. FALSE if it does not.
 */
BOOLEAN_T index_exists(
  INDEX_T* index ///< An allocated index
  )
{
  return index->on_disk;
}

/*
 * Private methods
 */

/*
 * Returns the index filename appropriate for this peptide
 */
char* get_peptide_file_name(
    INDEX_T* index,
    PEPTIDE_T* peptide
    );



/*********************************************
 * set and get methods for the object fields
 *********************************************/

/**
 *\returns the directory of the index
 * returns a heap allocated new copy of the directory
 * user must free the return directory name
 */
char* get_index_directory(
  INDEX_T* index ///< The index -in
  )
{
  return my_copy_string(index->directory);
}

/**
 * sets the directory of the index
 * index->directory must been initiailized
 */
void set_index_directory(
  INDEX_T* index, ///< The index -in
  const char* directory ///< the directory to add -in
  )
{
  free(index->directory);
  index->directory = my_copy_string(directory);
}

/**
 *\returns a pointer to the database
 */
DATABASE_T* get_index_database(
  INDEX_T* index ///< The index -in
  )
{
  return index->database;
}

/**
 * sets the database of the index
 */
void set_index_database(
  INDEX_T* index, ///< The index -in
  DATABASE_T* database ///< The database that has been indexed. -in
  )
{
  index->database = database;
}


/**
 * \brief Sets the peptide search constraint to be used by the
 * generate_peptides_iterator.  Makes a copy of the constraint pointer.
 * Deletes any existing search constraint. 
 */
void set_index_search_constraint(
  INDEX_T* index, ///< The index -in
  PEPTIDE_CONSTRAINT_T* constraint ///< Constraint for the next iterator
  )
{
  if(  index->search_constraint ){
    free_peptide_constraint(index->search_constraint);
  }
  index->search_constraint = copy_peptide_constraint_ptr(constraint);
  // check that the new mass window is within the index
  double search_min = get_peptide_constraint_min_mass(constraint);
  double search_max = get_peptide_constraint_max_mass(constraint);
  double index_min = get_peptide_constraint_min_mass(index->disk_constraint);
  double index_max = get_peptide_constraint_max_mass(index->disk_constraint);
  
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
 *\returns TRUE if index files are on disk else FALSE
 */
BOOLEAN_T get_index_on_disk(
  INDEX_T* index ///< The index -in
  )
{
  return index->on_disk;
}

/**
 * sets the on disk field of index
 */
void set_index_on_disk(
  INDEX_T* index, ///< The index -in
  BOOLEAN_T on_disk ///< Does this index exist on disk yet? -in
  )
{
  index->on_disk = on_disk;
}

/**
 *\returns the range of mass that each index file should be partitioned into
 */
FLOAT_T get_index_mass_range(
  INDEX_T* index ///< The index -in
  )
{
  return index->mass_range;
}

/**
 * sets the mass_range field of index
 */
void set_index_mass_range(
  INDEX_T* index, ///< The index -in
  FLOAT_T mass_range  ///< the range of mass that each index file should be partitioned into -in
  )
{
  index->mass_range = mass_range;
}

/**
 *\returns TRUE if only allow unique peptides else FALSE
 */
BOOLEAN_T get_index_is_unique(
  INDEX_T* index ///< The index -in
  )
{
  return index->is_unique;
}

/**
 * sets the is_unique field
 */
void set_index_is_unique(
  INDEX_T* index, ///< The index -in
  BOOLEAN_T is_unique ///< do you allow duplicate peptides? -in
  )
{
  index->is_unique = is_unique;
}


/**************************
 * Index file
 **************************/

// FIXME see if filename is a heap allocated
/**
 *\returns a new heap allocated index file object
 */
INDEX_FILE_T* new_index_file(
  char* filename,  ///< the filename to add -in
  FLOAT_T start_mass,  ///< the start mass of the index file  -in
  FLOAT_T range  ///< the mass range of the index file  -in
  )
{
  INDEX_FILE_T* index_file =
    (INDEX_FILE_T*)mycalloc(1, sizeof(INDEX_FILE_T));
  
  index_file->filename = filename;
  index_file->start_mass = start_mass;
  index_file->interval = range;

  return index_file;
}

/**
 * \brief Adds a new index_file object to the index_file.  Checks that
 * the total number of files does not exceed the limit.  Increases the
 * total_index_files count.
 * \returns TRUE if successfully added the new index_file
 */
BOOLEAN_T add_new_index_file(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator,
  ///< the index_peptide_iterator to add file -out
  char* filename_parsed,  ///< the filename to add -in
  FLOAT_T start_mass,  ///< the start mass of the index file  -in
  FLOAT_T range  ///< the mass range of the index file  -in
  )
{
  char* filename = my_copy_string(filename_parsed);
  carp(CARP_DETAILED_DEBUG, "Adding index file %s to iterator", filename);
  
  // check if total index files exceed MAX limit
  if(index_peptide_iterator->total_index_files > MAX_INDEX_FILES-1){
    carp(CARP_WARNING, "too many index files to read");
    return FALSE;
  }
  // create new index_file
  index_peptide_iterator->index_files[index_peptide_iterator->total_index_files] =
    new_index_file(filename, start_mass, range);
  
  ++index_peptide_iterator->total_index_files;
  return TRUE;
}

/**
 * frees the index file
 */
void free_index_file(
  INDEX_FILE_T* index_file  ///< the index file object to free -in
  )
{
  free(index_file->filename);
  free(index_file);
}

/*************************************
 * index_peptide iterator subroutines
 *************************************/

/**
 * \brief Parses the "crux_index_map" file that contains the mapping
 * between each crux_index_* file and a mass range. Adds all
 * crux_index_* files that are within the peptide constraint mass
 * range. 
 * \returns TRUE if successfully parses crux_index_map
 */
BOOLEAN_T parse_crux_index_map(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator  
    ///< working index_peptide_iterator -in
  )
{
  FILE* file = NULL;
  
  // used to parse each line from file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  
  // used to parse within a line

  FLOAT_T start_mass;
  FLOAT_T range;
  BOOLEAN_T start_file = FALSE;
  FLOAT_T min_mass = 
    get_peptide_constraint_min_mass(
          index_peptide_iterator->index->search_constraint);
  FLOAT_T max_mass = 
    get_peptide_constraint_max_mass(
          index_peptide_iterator->index->search_constraint);

  // used as buffer for reading in from file
  char full_filename[MAX_FILE_NAME_LENGTH] = "";
  strcpy(full_filename, index_peptide_iterator->index->directory);
  int dir_name_length = strlen(full_filename);

  // add a / to end of directory
  if( full_filename[dir_name_length-1] != '/' ){
    full_filename[dir_name_length] = '/';
    dir_name_length++;
  }
  // for filename as read from map file
  char* filename = full_filename + dir_name_length;
  // first use to open map file
  strcpy(filename, "crux_index_map");

  // open crux_index_file
  carp(CARP_DETAILED_DEBUG, "Opening map file '%s'", full_filename);
  file = fopen(full_filename, "r");
  if(file == NULL){
    int errsv = errno;
    carp(CARP_WARNING, "Cannot open crux_index_map file.:%s\nError:%s", 
      full_filename, 
      strerror(errsv));
    return FALSE;
  }
  
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    carp(CARP_DETAILED_DEBUG, "Index map file line reads %s", new_line);

    if(new_line[0] == 'c' && new_line[1] == 'r'){
      carp(CARP_DETAILED_DEBUG, "Looking for index file ");
      // read the crux_index_file information

      //      if(sscanf(new_line,"%s %f %f", 
      //                filename, &start_mass, &range) < 3){
      #ifdef USE_DOUBLES
      int char_read = sscanf(new_line,"%s %lf %lf", 
                             filename, &start_mass, &range);
      #else
      int char_read = sscanf(new_line,"%s %f %f", 
                             filename, &start_mass, &range);
      #endif
      if(char_read != 3){
        free(new_line);
        carp(CARP_WARNING, "Incorrect file format");
        fclose(file);
        return FALSE;
      }
      // find the first index file within mass range
      if(!start_file){
        if(min_mass > start_mass + range - 0.0001){
          continue;
        }
        else{
          start_file = TRUE;
          if(!add_new_index_file(
                index_peptide_iterator, full_filename, start_mass, range)){
            carp(CARP_WARNING, "Failed to add index file");
            fclose(file);
            free(new_line);
            return FALSE;
          }
          continue;
        }
      }// already added first file, add more
      // add all index_files that are with in peptide constraint mass interval
      else if(max_mass > (start_mass - 0.0001)){
        if(!add_new_index_file(
            index_peptide_iterator, full_filename, start_mass, range)){
          carp(CARP_WARNING, "Failed to add index file");
          free(new_line);
          return FALSE;
        }
        continue;
      }
      // out of mass range
      break;
    }
  }
  free(new_line);
  fclose(file);
  return TRUE;
}


/**
 * \brief Find the next peptide in the file that meets the peptide
 * constraint, read it from file, and add set up the iterator to
 * return it.
 * \returns TRUE if successfully finds and parses a peptide that meets
 * the constraint.
 */

BOOLEAN_T fast_forward_index_file(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator//, 
  ///< working index_peptide_iterator -in/out
  //FILE* file ///< the file stream to fast foward -in
  )
{
  FILE* file = index_peptide_iterator->index_file;
  // peptide to parse, reuse this memory while we look
  PEPTIDE_T* peptide = allocate_peptide();

  PEPTIDE_CONSTRAINT_T* index_constraint = 
    index_peptide_iterator->index->search_constraint;

  // loop until we get to a peptide that fits the constraint, we find
  // a peptide bigger (mass) than the constraint, or reach eof
  BOOLEAN_T peptide_fits = FALSE;
  long int src_loc = 0;
  while( ! peptide_fits ){
    // read in next peptide, returns false if eof
    if( ! parse_peptide_no_src(peptide, file, &src_loc) ){
      free_peptide(peptide);
      return FALSE;
    }
    // get mass & length
    int peptide_mass = (int)get_peptide_peptide_mass(peptide);
    int peptide_length = get_peptide_length(peptide);
    
    // if peptide mass larger than constraint, no more peptides to return
    if(peptide_mass > get_peptide_constraint_max_mass(index_constraint)){
      //free(peptide);
      free_peptide(peptide);
      return FALSE;
    }

    // does this peptide fit the constraint?
    double min_mass = get_peptide_constraint_min_mass(index_constraint);
    int min_length = get_peptide_constraint_min_length(index_constraint);
    int max_length = get_peptide_constraint_max_length(index_constraint);
    if( peptide_mass >= min_mass 
        && peptide_length >= min_length
        && peptide_length <= max_length){
      peptide_fits = TRUE;
    }
  } // read next peptide
  // now we have the peptide to hand to the iterator, finish parsing it

  // get peptide_src for this peptide
  long int pep_end = ftell(file);
  fseek(file, src_loc, SEEK_SET);
  DATABASE_T* database = index_peptide_iterator->index->database;
  if( ! parse_peptide_src(peptide, file, database, TRUE) ){
    carp(CARP_ERROR, "Could not parse peptide src");
    free_peptide(peptide);
    return FALSE;
  }

  fseek(file, pep_end, SEEK_SET);
  index_peptide_iterator->index_file = file;
  
  // add peptide to iterator
  index_peptide_iterator->peptide = peptide;

  return TRUE;
}

/**
 * clean_up
 * cleans up the temporary directory when SIGINT
 */
void clean_up( int dummy ) {

  fcloseall();
  delete_dir(temp_folder_name);
  exit(1);
  
  // quiet compiler
  dummy = dummy;
}


/**
 * \brief Prepare an index peptide iterator to have its index files
 * searched.  Changes the state of the iterator if its index_file
 * is NULL and the current file is not the last in the list.  The
 * routine that searches for peptides in a file closes any files it
 * has read to the end and sets the pointer to NULL.
 * \returns TRUE if there is a file ready to be read or FALSE if no
 * more files remain.
 */
BOOLEAN_T find_next_index_file(
  INDEX_PEPTIDE_ITERATOR_T* iterator
  ){
  carp(CARP_DETAILED_DEBUG, "Finding file");
  // file is ready to read
  if( iterator->index_file != NULL ){
    carp(CARP_DETAILED_DEBUG, "Current file ready to read");
    return TRUE;
  }
  // no more files to open
  if( iterator->current_index_file == iterator->total_index_files ){
  carp(CARP_DETAILED_DEBUG, "the last file has been opened. no more");
    return FALSE;
  }
  // else, current is NULL and there are more to open
  char* filename=iterator->index_files[iterator->current_index_file]->filename;
  carp(CARP_DETAILED_DEBUG, "Opening new file %s", filename);
  iterator->index_file = fopen(filename, "r");

  if( iterator->index_file == NULL){
    carp(CARP_ERROR, "Could not open index file %s", filename);
    return FALSE;
  }

  return TRUE;

}
  
/**
 * \brief Search for a peptide matching the peptide constraint in the
 * current index file of the iterator.  If the iterator has no open
 * index file, returns FALSE.  If the current file reaches the end
 * without finding a suitable peptide, the file is closed, file handle
 * set to NULL, and the number of the current file is increemnted.
 *
 * \returns TRUE if the iterator is ready to return a new peptide,
 * FALSE if no peptide meeting the constraint could abe found in the
 * current file. 
 */
BOOLEAN_T find_peptide_in_current_index_file(
  INDEX_PEPTIDE_ITERATOR_T* iterator)
{
  if( iterator == NULL){
    carp(CARP_ERROR, "Can't find peptide for NULL index iterator.");
    return FALSE;
  }
  carp(CARP_DETAILED_DEBUG, "Looking for peptide in current file");

  FILE* cur_file = iterator->index_file;
  if( cur_file == NULL ){
    carp(CARP_DETAILED_DEBUG, "current file is null");
    return FALSE;
  }

  // peptide to return, reuse this memory while we look
  PEPTIDE_T* peptide = allocate_peptide();
  // constraint to meet
  PEPTIDE_CONSTRAINT_T* index_constraint = iterator->index->search_constraint;

  // loop until we get to a peptide that fits the constraint, 
  // a peptide bigger (mass) than the constraint, or reach eof
  BOOLEAN_T peptide_fits = FALSE;
  BOOLEAN_T file_finished = FALSE;
  long int src_loc = 0;  // in case we need to parse the peptide src

  while( !peptide_fits && !file_finished ){// until pep_fits or file done
    carp(CARP_DETAILED_DEBUG, "Once around the find peptide loop %i, %i", 
         peptide_fits, file_finished);

    // read in next peptide
    BOOLEAN_T found_pep = parse_peptide_no_src(peptide, cur_file, &src_loc);
    // returns false if eof
    if( ! found_pep ){
      carp(CARP_DETAILED_DEBUG, "parse peptide returned FALSE");
      file_finished = TRUE;
      continue;
    }

    // check our peptide to see if it fits the constraint
    FLOAT_T peptide_mass = get_peptide_peptide_mass(peptide);
    int peptide_length = get_peptide_length(peptide);

    // if peptide mass larger than constraint, no more peptides to return
    if(peptide_mass > get_peptide_constraint_max_mass(index_constraint)){
      carp(CARP_DETAILED_DEBUG, "peptide found is bigger than constraint");
      file_finished = TRUE;
      continue;
    }// use else if?

    // does this peptide fit the constraint?

    double min_mass = get_peptide_constraint_min_mass(index_constraint);
    int min_length = get_peptide_constraint_min_length(index_constraint);
    int max_length = get_peptide_constraint_max_length(index_constraint);


    /*
  carp(CARP_DETAILED_DEBUG, "Index search constraints: %i-%i, %.2f-%.2f", get_peptide_constraint_min_length(iterator->index->search_constraint), get_peptide_constraint_max_length(iterator->index->search_constraint), get_peptide_constraint_min_mass(iterator->index->search_constraint), get_peptide_constraint_max_mass(iterator->index->search_constraint));
  carp(CARP_DETAILED_DEBUG, "Index search constraints var: %i-%i, %.2f-2f", min_length, max_length, min_mass); //, max_mass);
    */


    if( peptide_mass >= min_mass ){
      carp(CARP_DETAILED_DEBUG, "peptide mass bigger than min mass.");
    }else{
      carp(CARP_DETAILED_DEBUG, "peptide mass %f smaller than min mass %f.",
           peptide_mass, min_mass);
    }
    if( peptide_length >= min_length ){
      carp(CARP_DETAILED_DEBUG, "peptide length bigger than min length.");
    }
    if( peptide_length <= max_length ){
      carp(CARP_DETAILED_DEBUG, "peptide length bigger than min length.");
    }

    if( peptide_mass >= min_mass
        && peptide_length >= min_length
        && peptide_length <= max_length){
      carp(CARP_DETAILED_DEBUG, "peptide passes constraint.");
      peptide_fits = TRUE;
    }else{
      carp(CARP_DETAILED_DEBUG, "peptide doesn't pass constraint. " \
           "pep len,mass: %i, %.2f, min: %.2f, len: %i-%i", 
           peptide_length, peptide_mass, min_mass, min_length, max_length);
    }

  }// read next peptide

  // we have a peptide to return, get the peptide_src for it
  long int pep_end = 0;
  if( peptide_fits ){
    carp(CARP_DETAILED_DEBUG, "Found a peptide that fits constraint");
    pep_end = ftell(cur_file);
    fseek(cur_file, src_loc, SEEK_SET);
    DATABASE_T* database = iterator->index->database;
    if( ! parse_peptide_src(peptide, cur_file, database, TRUE) ){
      carp(CARP_ERROR, "Could not parse peptide src");
      file_finished = TRUE; // maybe we could read more, but unlikly
    }
  }

  // we broke out of the loop because we don't need to read the file anymore
  if( file_finished  ){
      carp(CARP_DETAILED_DEBUG, "Done with this index file.");
      fclose(cur_file);
      iterator->index_file = NULL;
      iterator->current_index_file += 1;
      iterator->has_next = FALSE;
      iterator->peptide = NULL;

      carp(CARP_DETAILED_DEBUG, "about to free peptide.");
      free_peptide(peptide);
      carp(CARP_DETAILED_DEBUG, "Done cleaning up.");
      return FALSE;
  }

  // else, everything worked!  finish setting fields

  //return file pointer to end of peptide
  fseek(cur_file, pep_end, SEEK_SET);
  iterator->index_file = cur_file;

  // add peptide to iterator
  iterator->peptide = peptide;
  iterator->has_next = TRUE;

  return TRUE;
}

/**
 * \brief Find the next peptide for the iterator to return.
 *
 * Called in the process of initializing a new iterator and by
 * next().  Looks in the current index_file for peptides matching the
 * constraint.  If none found, checks remaining index_files.  If no
 * constraint-satisfying peptides are found in any of the remaining
 * index files, next_peptide is set to NULL and has_next is set to FALSE.
 * \returns TRUE if no errors were encountered while reading files
 * (even if there is no peptide to return).
 */
BOOLEAN_T queue_next_peptide_index_peptide_iterator(
  INDEX_PEPTIDE_ITERATOR_T* iterator
  ){
  
  if(iterator == NULL){
    carp(CARP_ERROR, "Can't queue peptide for NULL index iterator.");
    return FALSE;
  }
  BOOLEAN_T found = FALSE;
  while(find_next_index_file(iterator)){
    found = find_peptide_in_current_index_file(iterator);
    if(found == TRUE ){
      carp(CARP_DETAILED_DEBUG, "Found returned TRUE");

      // set peptide, set has_next done in find
      return TRUE;
    }
  }// try the next index file

  // no more index files to try
  return FALSE;
}

/***********************************************
 *  The basic index_peptide_iterator functions.
 ***********************************************/

/**
 * Instantiates a new index_peptide_iterator from an index.
 * \returns a new heap allocated index_peptide_iterator object
 */
INDEX_PEPTIDE_ITERATOR_T* new_index_peptide_iterator(
  INDEX_T* index ///< The index object which we are iterating over -in
  )
{
  carp(CARP_DETAILED_DEBUG, "Creating new index iterator");

  // set peptide implementation to array peptide_src
  // this determines which peptide free method to use
  set_peptide_src_implementation(FALSE);

  // allocate a new index_peptide_iterator object
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator =
    (INDEX_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(INDEX_PEPTIDE_ITERATOR_T));
  
  // set index
  index_peptide_iterator->index = copy_index_ptr(index);
  
  // parse index_files that are within peptide_constraint from crux_index_map
  // sets index_files and total_index_files
  int parse_count = 0;
  while(!parse_crux_index_map(index_peptide_iterator)){
    // failed to parse crux_index_map
    if (parse_count++ > MAX_PARSE_COUNT){
      carp(CARP_FATAL, 
        "Failed to parse crux_index_map file after %i tries", MAX_PARSE_COUNT);
    } else {
      carp(CARP_ERROR, "Failed to parse crux_index_map file. Sleeping.");
      sleep(SLEEP_DURATION);
    }
  }

  // set remaining iterator fields
  index_peptide_iterator->current_index_file = 0;
  index_peptide_iterator->index_file = NULL;

  // sets has_next, peptide, index_file
  carp(CARP_DETAILED_DEBUG, "Queueing first peptide");
  queue_next_peptide_index_peptide_iterator(index_peptide_iterator);

  return index_peptide_iterator;
}

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* index_peptide_iterator_next(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator ///< the index_peptide_iterator to initialize -in
  )
{
  PEPTIDE_T* peptide_to_return = index_peptide_iterator->peptide;

  queue_next_peptide_index_peptide_iterator(index_peptide_iterator);

  return peptide_to_return;
}

/**
 * The basic iterator functions.
 * check to see if the index_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T index_peptide_iterator_has_next(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator ///< the index_peptide_iterator to initialize -in
  )
{
  return index_peptide_iterator->has_next; 
}

/**
 * Frees an allocated index_peptide_iterator object.
 */
void free_index_peptide_iterator(
  INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator ///< the iterator to free -in
  )
{
  
  // free all index files
  int file_idx;
  for(file_idx=0;file_idx<index_peptide_iterator->total_index_files;++file_idx){
    free_index_file(index_peptide_iterator->index_files[file_idx]);
  }
  
  // if did not iterate over all peptides, free the last peptide not returned
  if(index_peptide_iterator_has_next(index_peptide_iterator)){
    free_peptide(index_peptide_iterator->peptide);
  }
  
  // free the index
  free_index(index_peptide_iterator->index);
  
  free(index_peptide_iterator);
}

/***********************************************
 * index_filtered_peptide_iterator
 ***********************************************/

/**
 * sets up the index_filtered_peptide_iterator
 * \returns TRUE if successfully sets up the iterator, else FALSE
 */
BOOLEAN_T setup_index_filtered_peptide_iterator(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* iterator
  )
{
  PEPTIDE_T* peptide = NULL;
  PEPTIDE_SRC_T* src = NULL;
  /*
  PEPTIDE_TYPE_T peptide_type = 
    get_peptide_constraint_peptide_type(iterator->index_peptide_iterator->index->search_constraint);
  */
  DIGEST_T required_digestion = get_peptide_constraint_digest(iterator->index_peptide_iterator->index->search_constraint);
  BOOLEAN_T match = FALSE;

  // initialize index_filered
  while(index_peptide_iterator_has_next(iterator->index_peptide_iterator)){
    peptide = index_peptide_iterator_next(iterator->index_peptide_iterator);
    src = get_peptide_peptide_src(peptide);
    // mass, length has been already checked in index_peptide_iterator
    // check if peptide type matches the constraint
    // find at least one peptide_src for which cleavage is correct
    while(src != NULL){
      /*
      char this_type_str[16];
      char request_type_str[16];
      peptide_type_to_string(peptide_type, request_type_str);
      peptide_type_to_string(get_peptide_src_peptide_type(src), this_type_str);
      carp(CARP_DETAILED_DEBUG, "request: %s, this: %s", request_type_str, this_type_str);
      if(get_peptide_src_peptide_type(src) == peptide_type ||
         (peptide_type == PARTIALLY_TRYPTIC && 
          (get_peptide_src_peptide_type(src) == N_TRYPTIC ||
           get_peptide_src_peptide_type(src) == C_TRYPTIC))
         ){
      */
      if(get_peptide_src_digest(src) >= required_digestion){
        match = TRUE;
        break;
      }
      // check the next peptide src
      src = get_peptide_src_next_association(src);
    }

    // add more filters to the peptides here, if they don't meet
    // requirements change 'match' to FALSE 

    // this peptide meets the peptide_type
    if(match){
      iterator->peptide = peptide;
      iterator->has_next = TRUE;
      return TRUE;
    }
    free_peptide(peptide);

  }// next peptide
  // no peptides meet the constraint
  iterator->has_next = FALSE;
  return TRUE;
}

/**
 * Instantiates a new index_filtered_peptide_iterator from a index.
 * \returns a new heap allocated index_filtered_peptide_iterator object
 */
INDEX_FILTERED_PEPTIDE_ITERATOR_T* new_index_filtered_peptide_iterator(
  INDEX_T* index ///< The index object which we are iterating over -in
  )
{
  // create new index_filtered_peptide_iterator
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_iterator =
    (INDEX_FILTERED_PEPTIDE_ITERATOR_T*)mycalloc(
        1, sizeof(INDEX_FILTERED_PEPTIDE_ITERATOR_T));

  // create new index peptide iterator, the core peptide iterator
  index_filtered_iterator->index_peptide_iterator = 
    new_index_peptide_iterator(index);
  
  // setup index_filered_iterator
  if(!setup_index_filtered_peptide_iterator(index_filtered_iterator)){
    carp(CARP_ERROR, "Failed to setup index filtered peptide iterator");
  }
  
  return index_filtered_iterator;
}

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* index_filtered_peptide_iterator_next(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  PEPTIDE_T* peptide_to_return = index_filtered_peptide_iterator->peptide;

  // check if there's actually a peptide to return
  if(!index_filtered_peptide_iterator_has_next(index_filtered_peptide_iterator)
     || index_filtered_peptide_iterator->peptide == NULL){
    return NULL;
  }
  
  // setup the interator for the next peptide, if avaliable
  if(!setup_index_filtered_peptide_iterator(index_filtered_peptide_iterator)){
    carp(CARP_FATAL, "failed to setup index_filtered_peptide_iterator for next iteration");
  }
  return peptide_to_return;
}

/**
 * The basic iterator functions.
 * check to see if the index_filtered_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T index_filtered_peptide_iterator_has_next(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  return index_filtered_peptide_iterator->has_next; 
}

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void free_index_filtered_peptide_iterator(
  INDEX_FILTERED_PEPTIDE_ITERATOR_T* index_filtered_peptide_iterator ///< the iterator to free -in
  )
{
  free_index_peptide_iterator(
      index_filtered_peptide_iterator->index_peptide_iterator);
    
  // if did not iterate over all peptides, free the last peptide not returned
  if(index_filtered_peptide_iterator_has_next(index_filtered_peptide_iterator)){
    free_peptide(index_filtered_peptide_iterator->peptide);
  }

  free(index_filtered_peptide_iterator);
}

/****************************
 * bin_peptide_iterator
 ******************************/

/**
 * reads in the first peptide from the bin
 * parse the peptide, then adds it to the bin-peptide-iterator to return
 * \returns TRUE if successfully initializes the bin_peptide_iterator
 */
BOOLEAN_T initialize_bin_peptide_iterator(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< working bin_peptide_iterator -in/out
  )
{
  // BUG: This no longer distinguishes between eof and error in
  // parsing.  One fix would be for parse_peptide to return error code

  FILE* file = bin_peptide_iterator->index_file;
  DATABASE_T* database = bin_peptide_iterator->index->database;
  BOOLEAN_T use_src_array = bin_peptide_iterator->use_array;

  // allocate peptide to used to parse
  //  PEPTIDE_T* peptide = allocate_peptide();
  PEPTIDE_T* peptide = parse_peptide(file, database, use_src_array);

  if( peptide == NULL ){
    bin_peptide_iterator->peptide = NULL;
    bin_peptide_iterator->has_next = FALSE;
    //return FALSE;
    return TRUE;
  }
  // set file pointer
  bin_peptide_iterator->index_file = file;
  bin_peptide_iterator->peptide = peptide;
  bin_peptide_iterator->has_next = TRUE;
  return TRUE;
}


/**
 * Instantiates a new bin_peptide_iterator from a gvien bin file handler.
 * \returns a new heap allocated bin_peptide_iterator object
 */
BIN_PEPTIDE_ITERATOR_T* new_bin_peptide_iterator(
  INDEX_T* index, ///< The index object which we are iterating over -in
  FILE* file, ///< the bin to parse peptides
  BOOLEAN_T use_array  ///< should I use array peptide_src or link list when parsing peptides -in
  )
{
  if(use_array){
    // set peptide implementation to array peptide_src
    // this determines which peptide free method to use
    set_peptide_src_implementation(FALSE);
  }
  else{// use link list
    set_peptide_src_implementation(TRUE);
  }

  // allocate a new index_peptide_iterator object
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator =
    (BIN_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(BIN_PEPTIDE_ITERATOR_T));
  
  // set index, file
  bin_peptide_iterator->index = index;
  bin_peptide_iterator->index_file = file;
  bin_peptide_iterator->use_array = use_array;
    
  if(!initialize_bin_peptide_iterator(bin_peptide_iterator)){
    carp(CARP_WARNING, "failed to initalize bin peptide iterator");
    bin_peptide_iterator->has_next = FALSE;
  }
    
  return bin_peptide_iterator;
}


/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* bin_peptide_iterator_next(
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator ///< the bin_peptide_iterator to get peptide -in
  )
{
  PEPTIDE_T* peptide_to_return = bin_peptide_iterator->peptide;

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
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T bin_peptide_iterator_has_next(
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
    free_peptide(bin_peptide_iterator->peptide);
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
  INDEX_T* index, ///< The index object which we are iterating over -in
  FILE* file, ///< the working bin file handler -in
  unsigned int peptide_count ///< the total peptide count in the bin -in
  )
{
  // set peptide implementation to array peptide_src
  // this determines which peptide free method to use
  set_peptide_src_implementation(FALSE);
  
  // create database sorted peptide iterator
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator =
    (BIN_SORTED_PEPTIDE_ITERATOR_T*)
    mycalloc(1, sizeof(BIN_SORTED_PEPTIDE_ITERATOR_T));

  // reset file to start
  rewind(file);

  // create bin_peptide_iterator
  // use link list peptide_src implementation to merge peptides
  BIN_PEPTIDE_ITERATOR_T* bin_peptide_iterator =
    new_bin_peptide_iterator(index, file, FALSE);

  // create a sorted peptide iterator that will sort all the peptides 
  // from bin peptide_iterator
  SORTED_PEPTIDE_ITERATOR_T* sorted_peptide_iterator = 
    new_sorted_peptide_iterator_bin(bin_peptide_iterator, SORT_MASS, 
        index->is_unique, peptide_count);

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
PEPTIDE_T* bin_sorted_peptide_iterator_next(
  BIN_SORTED_PEPTIDE_ITERATOR_T* bin_sorted_peptide_iterator 
    ///< the bin_peptide_iterator to get peptide -in
  )
{
  PEPTIDE_T* peptide = sorted_peptide_iterator_next(bin_sorted_peptide_iterator->sorted_peptide_iterator);
  return peptide;
}

/**
 * The basic iterator functions.
 * check to see if the bin_sorted_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T bin_sorted_peptide_iterator_has_next(
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



/**********************************************************************
 * wrapper, for generate_peptides_iterator, cast back to original type
 ***********************************************************************/

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* void_index_filtered_peptide_iterator_next(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  return index_filtered_peptide_iterator_next((INDEX_FILTERED_PEPTIDE_ITERATOR_T*)index_filtered_peptide_iterator);

}

/**
 * The basic iterator functions.
 * check to see if the index_filtered_peptide_iterator has more peptides to return
 *\returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_index_filtered_peptide_iterator_has_next(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  return index_filtered_peptide_iterator_has_next((INDEX_FILTERED_PEPTIDE_ITERATOR_T*)index_filtered_peptide_iterator);
}

/**
 * Frees an allocated index_filtered_peptide_iterator object.
 */
void void_free_index_filtered_peptide_iterator(
  void* index_filtered_peptide_iterator ///< the index_filtered_peptide_iterator to initialize -in
  )
{
  free_index_filtered_peptide_iterator(
      (INDEX_FILTERED_PEPTIDE_ITERATOR_T*)index_filtered_peptide_iterator);
}

/**
 * Frees an allocated index_peptide_iterator object.
 */
void void_free_index_peptide_iterator(
    void* index_peptide_iterator  ///< the iterator to free -in
    )
{
  free_index_peptide_iterator(
      (INDEX_PEPTIDE_ITERATOR_T*)index_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T void_index_peptide_iterator_has_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    )
{

  return index_peptide_iterator_has_next((INDEX_PEPTIDE_ITERATOR_T*)index_peptide_iterator);
}

/**
 * \returns The next peptide in the index.
 */
PEPTIDE_T* void_index_peptide_iterator_next(
    void* index_peptide_iterator ///< the iterator of interest -in
    )
{

  return index_peptide_iterator_next((INDEX_PEPTIDE_ITERATOR_T*)index_peptide_iterator);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
