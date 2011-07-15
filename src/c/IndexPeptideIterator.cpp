/************************************************************************//**
 * \file IndexPeptideIterator.cpp
 * \brief Object for representing iterating peptides from an index.
 ****************************************************************************/

#include "IndexPeptideIterator.h"
#include <errno.h>

static const int MAX_FILE_NAME_LENGTH = 300;
static const int MAX_PARSE_COUNT = 3;
static const int SLEEP_DURATION = 5;


/***********************************************
 *  The basic index_peptide_iterator functions.
 ***********************************************/

/**
 * Instantiates a new index_peptide_iterator from an index.
 * \returns a new heap allocated index_peptide_iterator object
 */
IndexPeptideIterator::IndexPeptideIterator(
  Index* index ///< The index object which we are iterating over -in
  )
{
  carp(CARP_DETAILED_DEBUG, "Creating new index iterator");

  // set peptide implementation to array peptide_src
  // this determines which peptide free method to use
  set_peptide_src_implementation(false);

  // initialize a new index_peptide_iterator object
  index_ = NULL;
  total_index_files_ = 0;
  current_index_file_ = 0;
  index_file_ = NULL;
  has_next_ = false;
  peptide_ = NULL;
  
  // set index
  index_ = index->copyPtr();
  
  // parse index_files that are within peptide_constraint from crux_index_map
  // sets index_files and total_index_files
  int parse_count = 0;
  while(!parseCruxIndexMap()){
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
  current_index_file_ = 0;
  index_file_ = NULL;

  // sets has_next, peptide, index_file
  carp(CARP_DETAILED_DEBUG, "Queueing first peptide");
  queueNextPeptide();
}

/**
 *  The basic iterator functions.
 * \returns The next peptide in the index.
 */
PEPTIDE_T* IndexPeptideIterator::next()
{
  PEPTIDE_T* peptide_to_return = peptide_;

  queueNextPeptide();

  return peptide_to_return;
}

/**
 * The basic iterator functions.
 * check to see if the index_peptide_iterator has more peptides to return
 *\returns true if there are additional peptides to iterate over, false if not.
 */
bool IndexPeptideIterator::hasNext()
{
  return has_next_; 
}

/**
 * Frees an allocated index_peptide_iterator object.
 */
IndexPeptideIterator::~IndexPeptideIterator()
{
  
  // free all index files
  int file_idx;
  for(file_idx=0;file_idx<total_index_files_;++file_idx){
    delete index_files_[file_idx];
  }
  
  // if did not iterate over all peptides, free the last peptide not returned
  if(hasNext()){
    free_peptide(peptide_);
  }
  
  // free the index
  Index::free(index_);
}

/**
 * \brief Prepare an index peptide iterator to have its index files
 * searched.  Changes the state of the iterator if its index_file
 * is NULL and the current file is not the last in the list.  The
 * routine that searches for peptides in a file closes any files it
 * has read to the end and sets the pointer to NULL.
 * \returns true if there is a file ready to be read or false if no
 * more files remain.
 */
bool IndexPeptideIterator::findNextIndexFile(){

  carp(CARP_DETAILED_DEBUG, "Finding file");
  // file is ready to read
  if( index_file_ != NULL ){
    carp(CARP_DETAILED_DEBUG, "Current file ready to read");
    return true;
  }
  // no more files to open
  if( current_index_file_ == total_index_files_ ){
  carp(CARP_DETAILED_DEBUG, "the last file has been opened. no more");
    return false;
  }
  // else, current is NULL and there are more to open
  char* filename=index_files_[current_index_file_]->getFilename();
  carp(CARP_DETAILED_DEBUG, "Opening new file %s", filename);
  index_file_ = fopen(filename, "r");

  if( index_file_ == NULL){
    carp(CARP_ERROR, "Could not open index file %s", filename);
    return false;
  }

  return true;

}
  


/**
 * \brief Search for a peptide matching the peptide constraint in the
 * current index file of the iterator.  If the iterator has no open
 * index file, returns false.  If the current file reaches the end
 * without finding a suitable peptide, the file is closed, file handle
 * set to NULL, and the number of the current file is increemnted.
 *
 * \returns true if the iterator is ready to return a new peptide,
 * false if no peptide meeting the constraint could abe found in the
 * current file. 
 */
bool IndexPeptideIterator::findPeptideInCurrentIndexFile()
{
  carp(CARP_DETAILED_DEBUG, "Looking for peptide in current file");

  FILE* cur_file = index_file_;
  if( cur_file == NULL ){
    carp(CARP_DETAILED_DEBUG, "current file is null");
    return false;
  }

  // peptide to return, reuse this memory while we look
  PEPTIDE_T* peptide = allocate_peptide();
  // constraint to meet
  PeptideConstraint* index_constraint = index_->getSearchConstraint();

  // loop until we get to a peptide that fits the constraint, 
  // a peptide bigger (mass) than the constraint, or reach eof
  bool peptide_fits = false;
  bool file_finished = false;
  long int src_loc = 0;  // in case we need to parse the peptide src

  while( !peptide_fits && !file_finished ){// until pep_fits or file done
    carp(CARP_DETAILED_DEBUG, "Once around the find peptide loop %i, %i", 
         peptide_fits, file_finished);

    // read in next peptide
    bool found_pep = parse_peptide_no_src(peptide, cur_file, &src_loc);
    // returns false if eof
    if( ! found_pep ){
      carp(CARP_DETAILED_DEBUG, "parse peptide returned false");
      file_finished = true;
      continue;
    }

    // check our peptide to see if it fits the constraint
    FLOAT_T peptide_mass = get_peptide_peptide_mass(peptide);
    int peptide_length = get_peptide_length(peptide);

    // if peptide mass larger than constraint, no more peptides to return
    if(peptide_mass > index_constraint->getMaxMass()){
      carp(CARP_DETAILED_DEBUG, "peptide found is bigger than constraint");
      file_finished = true;
      continue;
    }// use else if?

    // does this peptide fit the constraint?

    double min_mass = index_constraint->getMinMass();
    int min_length = index_constraint->getMinLength();
    int max_length = index_constraint->getMaxLength();


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
      carp(CARP_DETAILED_DEBUG, "peptide length smaller than max length.");
    }

    if( peptide_mass >= min_mass
        && peptide_length >= min_length
        && peptide_length <= max_length){
      carp(CARP_DETAILED_DEBUG, "peptide passes constraint.");
      peptide_fits = true;
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
    Database* database = index_->getDatabase();
    if( ! parse_peptide_src(peptide, cur_file, database, true) ){
      carp(CARP_ERROR, "Could not parse peptide src");
      file_finished = true; // maybe we could read more, but unlikly
    }
  }

  // we broke out of the loop because we don't need to read the file anymore
  if( file_finished  ){
      carp(CARP_DETAILED_DEBUG, "Done with this index file.");
      fclose(cur_file);
      index_file_ = NULL;
      current_index_file_ += 1;
      has_next_ = false;
      peptide_ = NULL;

      carp(CARP_DETAILED_DEBUG, "about to free peptide.");
      free_peptide(peptide);
      carp(CARP_DETAILED_DEBUG, "Done cleaning up.");
      return false;
  }

  // else, everything worked!  finish setting fields

  //return file pointer to end of peptide
  fseek(cur_file, pep_end, SEEK_SET);
  index_file_ = cur_file;

  // add peptide to iterator
  peptide_ = peptide;
  has_next_ = true;

  return true;
}

/**
 * \brief Find the next peptide for the iterator to return.
 *
 * Called in the process of initializing a new iterator and by
 * next().  Looks in the current index_file for peptides matching the
 * constraint.  If none found, checks remaining index_files.  If no
 * constraint-satisfying peptides are found in any of the remaining
 * index files, next_peptide is set to NULL and has_next is set to false.
 * \returns true if no errors were encountered while reading files
 * (even if there is no peptide to return).
 */
bool IndexPeptideIterator::queueNextPeptide() {

  bool found = false;
  while(findNextIndexFile()){
    found = findPeptideInCurrentIndexFile();
    if(found == true ){
      carp(CARP_DETAILED_DEBUG, "Found returned true");

      // set peptide, set has_next done in find
      return true;
    }
  }// try the next index file

  // no more index files to try
  return false;
}



/**
 * \brief Adds a new index_file object to the index_file.  Checks that
 * the total number of files does not exceed the limit.  Increases the
 * total_index_files count.
 * \returns true if successfully added the new index_file
 */
bool IndexPeptideIterator::addNewIndexFile(
  ///< the index_peptide_iterator to add file -out
  char* filename_parsed,  ///< the filename to add -in
  FLOAT_T start_mass,  ///< the start mass of the index file  -in
  FLOAT_T range  ///< the mass range of the index file  -in
  )
{
  char* filename = my_copy_string(filename_parsed);
  carp(CARP_DETAILED_DEBUG, "Adding index file %s to iterator", filename);
  
  // check if total index files exceed MAX limit
  if(total_index_files_ > MAX_INDEX_FILES-1){
    carp(CARP_WARNING, "too many index files to read");
    return false;
  }
  // create new index_file
  index_files_[total_index_files_] =
    new IndexFile(filename, start_mass, range);
  
  ++total_index_files_;
  return true;
}

/**
 * \brief Parses the "crux_index_map" file that contains the mapping
 * between each crux_index_* file and a mass range. Adds all
 * crux_index_* files that are within the peptide constraint mass
 * range. 
 * \returns true if successfully parses crux_index_map
 */
bool IndexPeptideIterator::parseCruxIndexMap()
{
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
    index_->getSearchConstraint()->getMinMass();
  FLOAT_T max_mass = 
    index_->getSearchConstraint()->getMaxMass();

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
  strcpy(filename, "crux_index_map");

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
 * \brief Find the next peptide in the file that meets the peptide
 * constraint, read it from file, and add set up the iterator to
 * return it.
 * \returns true if successfully finds and parses a peptide that meets
 * the constraint.
 */

bool IndexPeptideIterator::fastForwardIndexFile(
  ///< working index_peptide_iterator -in/out
  //FILE* file ///< the file stream to fast foward -in
  )
{
  FILE* file = index_file_;
  // peptide to parse, reuse this memory while we look
  PEPTIDE_T* peptide = allocate_peptide();

  PeptideConstraint* index_constraint = 
    index_->getSearchConstraint();

  // loop until we get to a peptide that fits the constraint, we find
  // a peptide bigger (mass) than the constraint, or reach eof
  bool peptide_fits = false;
  long int src_loc = 0;
  while( ! peptide_fits ){
    // read in next peptide, returns false if eof
    if( ! parse_peptide_no_src(peptide, file, &src_loc) ){
      free_peptide(peptide);
      return false;
    }
    // get mass & length
    int peptide_mass = (int)get_peptide_peptide_mass(peptide);
    int peptide_length = get_peptide_length(peptide);
    
    // if peptide mass larger than constraint, no more peptides to return
    if(peptide_mass > index_constraint->getMaxMass()){
      //free(peptide);
      free_peptide(peptide);
      return false;
    }

    // does this peptide fit the constraint?
    double min_mass = index_constraint->getMinMass();
    int min_length = index_constraint->getMinLength();
    int max_length = index_constraint->getMaxLength();
    if( peptide_mass >= min_mass 
        && peptide_length >= min_length
        && peptide_length <= max_length){
      peptide_fits = true;
    }
  } // read next peptide
  // now we have the peptide to hand to the iterator, finish parsing it

  // get peptide_src for this peptide
  long int pep_end = ftell(file);
  fseek(file, src_loc, SEEK_SET);
  Database* database = index_->getDatabase();
  if( ! parse_peptide_src(peptide, file, database, true) ){
    carp(CARP_ERROR, "Could not parse peptide src");
    free_peptide(peptide);
    return false;
  }

  fseek(file, pep_end, SEEK_SET);
  index_file_ = file;
  
  // add peptide to iterator
  peptide_ = peptide;

  return true;
}

Index* IndexPeptideIterator::getIndex() {
  
  return index_;
}
