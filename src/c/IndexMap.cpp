/************************************************************************//**
 * \file IndexMap.cpp
 * \brief Object for representing an the IndexMap of a database
 ****************************************************************************/

#include "IndexMap.h"
#include <errno.h>
#include "WinCrux.h"

#include "Index.h"

using namespace std;

static const int MAX_FILE_NAME_LENGTH = 300;
static const int MAX_PARSE_COUNT = 3;
static const int SLEEP_DURATION = 5;

/**
 * Loads the index map into memory
 */
IndexMap::IndexMap(
  Index* index, ///< Index where the map should be generated from
  bool is_decoy ///< Are we using the decoy or target index?
) {

  index_ = index;
  is_decoy_ = is_decoy;

  int parse_count = 0;
  

  while(!parse()){
    // failed to parse crux_index_map
    if (parse_count++ > MAX_PARSE_COUNT){
      carp(CARP_FATAL, 
        "Failed to parse crux_index_map file after %i tries", MAX_PARSE_COUNT);
    } else {
      carp(CARP_ERROR, "Failed to parse crux_index_map file. Sleeping.");
      sleep(SLEEP_DURATION);
    }
  }
}

/**
 * Default destructor
 */
IndexMap::~IndexMap() {

  for (size_t idx = 0; idx < index_files_.size() ; idx++) {
    delete index_files_[idx];
  }
  index_files_.clear();

}


/**
 * \brief Parses the "crux_index_map" file that contains the mapping
 * between each crux_index_* file and a mass range. Adds all
 * crux_index_* files that are within the disk constraint mass
 * range. 
 * \returns true if successfully parses crux_index_map
 */
bool IndexMap::parse() {
  FILE* file = NULL;
  
  // used to parse each line from file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  
  // used to parse within a line

  FLOAT_T start_mass;
  FLOAT_T range;
  bool start_file = false;
  FLOAT_T min_mass = 
    index_->getDiskConstraint()->getMinMass();
  FLOAT_T max_mass = 
    index_->getDiskConstraint()->getMaxMass();

  // used as buffer for reading in from file
  char full_filename[MAX_FILE_NAME_LENGTH] = "";
  strcpy(full_filename, index_->getDirectory());
  int dir_name_length = strlen(full_filename);

  // add a / to end of directory
  if( full_filename[dir_name_length-1] != '/' ){
    full_filename[dir_name_length] = '/';
    dir_name_length++;
  }
  // for filename as read from map file
  char* filename = full_filename + dir_name_length;
  // first use to open map file
  if( is_decoy_ ){
    strcpy(filename, Index::decoy_index_file_prefix);
    strcpy(filename + strlen(Index::decoy_index_file_prefix), "map");
  } else {
    strcpy(filename, Index::index_file_prefix);
    strcpy(filename + strlen(Index::index_file_prefix), "map");
  }

  // open crux_index_file
  carp(CARP_DETAILED_DEBUG, "Opening map file '%s'", full_filename);
  file = fopen(full_filename, "r");
  if(file == NULL){
    int errsv = errno;
    carp(CARP_WARNING, "Cannot open crux_index_map file.:%s\nError:%s", 
      full_filename, 
      strerror(errsv));
    return false;
  }
  
  while((line_length =  getline(&new_line, &buf_length, file)) != -1){
    carp(CARP_DETAILED_DEBUG, "Index map file line reads '%s'", new_line);

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
        return false;
      }
      // find the first index file within mass range
      if(!start_file){
        if(min_mass > start_mass + range - 0.0001){
          continue;
        }
        else{
          start_file = true;
          if(!addNewIndexFile(full_filename, start_mass, range)){
            carp(CARP_WARNING, "Failed to add index file");
            fclose(file);
            free(new_line);
            return false;
          }
          continue;
        }
      }// already added first file, add more
      // add all index_files that are with in peptide constraint mass interval
      else if(max_mass > (start_mass - 0.0001)){
        if(!addNewIndexFile(
            full_filename, start_mass, range)){
          carp(CARP_WARNING, "Failed to add index file");
          free(new_line);
          return false;
        }
        continue;
      }
      // out of mass range
      break;
    }
  }
  free(new_line);
  fclose(file);
  return true;
}

/**
 * \brief Adds a new index_file object to the index_file.  
 * \returns true if successfully added the new index_file
 */
bool IndexMap::addNewIndexFile(
  const char* filename_parsed,  ///< the filename to add -in
  FLOAT_T start_mass,  ///< the start mass of the index file  -in
  FLOAT_T range  ///< the mass range of the index file  -in
  )
{
  char* filename = my_copy_string(filename_parsed);
  carp(CARP_DETAILED_DEBUG, "Adding index file %s to iterator", filename);

  // create new index_file
  index_files_.push_back(new IndexFile(filename, start_mass, range));
  return true;
}

/**
 * populates the constraint_index_files with the index files that match
 * the passed in constraint
 */
void IndexMap::getIndexFiles(
  PeptideConstraint* constraint, ///< the constraint for peptides to select -in 
  vector<IndexFile*>& constraint_index_files ///< a list of all of the indices passing the constraint -out
  ) {

  constraint_index_files.clear();

  FLOAT_T min_mass = constraint->getMinMass();
  FLOAT_T max_mass = constraint->getMaxMass();

  for (size_t idx = 0; idx < index_files_.size();idx++) {
    IndexFile* index_file = index_files_[idx];  
    FLOAT_T start_mass = index_file->getStartMass();
    FLOAT_T range = index_file->getRange();

    if ((start_mass - 0.0001) >= max_mass) {
      break;
    }

    if ((start_mass + range - 0.0001) > min_mass) {
      constraint_index_files.push_back(index_file);
    }
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
