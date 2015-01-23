#include "crux-file-utils.h"


/**
 * \returns True if there already exists a file with the given name,
 * else false.
 */
bool file_exists(const char* filename){
  struct stat file_info;
  int stat_return = stat(filename, &file_info);
  if( stat_return == 0 ){
    return true;
  } // else, we could not get attributes.  It doesn't exist or we
    // don't have permission for the directory
  return false;
}

/**
 * Open a file of the given name if it either does not exist or if we
 * have permission to overwrite.
 * \returns A file stream for the open file or NULL on failure.
 */
std::ofstream* create_file
(const char* filename, ///< create file with this name
 bool overwrite ///< replace any existing files with this name.
 ){

  if( file_exists(filename) && !overwrite ){
    // print warning?
    return NULL;
  }

  std::ofstream* file = new std::ofstream(filename, std::ios::out);
  if( ! file->is_open() ){
    delete file;
    return NULL;
  }

  return file;
}


