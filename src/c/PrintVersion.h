/**
 * \file PrintVersion.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#ifndef PRINTVERSION_H
#define PRINTVERSION_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class PrintVersion: public CruxApplication {

 public:
  /**
   * \returns a blank PrintVersion object
   */
  PrintVersion();
  
  /**
   * Destructor
   */
  ~PrintVersion();

  /**
   * main method for PrintVersion
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for PrintVersion
   */
  virtual std::string getName();

  /**
   * \returns the description for PrintVersion
   */
  virtual std::string getDescription();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
