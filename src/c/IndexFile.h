#ifndef INDEXFILE_H
#define INDEXFILE_H

#include "utils.h"

class IndexFile {
  
 protected:
  char* filename_;  ///< The file name that contain the peptides
  FLOAT_T start_mass_; ///< the start mass limit in this file
  FLOAT_T interval_;   ///< the interval of the peptides in this file
 public:
  
  /**
   *\returns a new heap allocated index file object
   */
  IndexFile(
    char* filename,  ///< the filename to add -in
    FLOAT_T start_mass,  ///< the start mass of the index file  -in
    FLOAT_T range  ///< the mass range of the index file  -in
    );

  /**
   * frees the index file
   */
  virtual ~IndexFile();

  
  char* getFilename();
  FLOAT_T getStartMass();
  FLOAT_T getRange();

};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
