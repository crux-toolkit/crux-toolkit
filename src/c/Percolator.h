/**
 * \file Percolator.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running percolator
 *****************************************************************************/
#ifndef PERCOLATOR_H
#define PERCOLATOR_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class Percolator: public CruxApplication {

 public:
  /**
   * \returns a blank Percolator object
   */
  Percolator();
  
  /**
   * Destructor
   */
  ~Percolator();

  /**
   * main method for Percolator
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for Percolator
   */
  virtual std::string getName();

  /**
   * \returns the description for Percolator
   */
  virtual std::string getDescription();
  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem();

  virtual COMMAND_T getCommand();
  
  /**
   * \returns whether the application needs the output directory or not.
   */
  virtual bool needsOutputDirectory();


};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

