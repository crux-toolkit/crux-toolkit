#include "IndexFile.h"


/**************************
 * Index file
 **************************/

// FIXME see if filename is a heap allocated
/**
 *\returns a new heap allocated index file object
 */
IndexFile::IndexFile(
  char* filename,  ///< the filename to add -in
  FLOAT_T start_mass,  ///< the start mass of the index file  -in
  FLOAT_T range  ///< the mass range of the index file  -in
  )
{

  filename_ = filename;
  start_mass_ = start_mass;
  interval_ = range;
}


/**
 * frees the index file
 */
IndexFile::~IndexFile() {

  free(filename_);

}


char* IndexFile::getFilename() {
  return filename_;
}

FLOAT_T IndexFile::getStartMass() {
  return start_mass_;

}

FLOAT_T IndexFile::getRange() {
  return interval_;
}
