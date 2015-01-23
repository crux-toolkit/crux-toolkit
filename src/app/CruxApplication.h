/**
 * \file CruxApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#ifndef CRUXAPPLICATION_H
#define CRUXAPPLICATION_H
#include "objects.h"


#include <string>

class CruxApplication{
 public:
  /**
   * the main method for CruxApplication.  Subclasses of
   * CruxApplication define this.
   * \returns exit code for the executed program.   
   */
  virtual int main(int argc, char** argv)=0;

  /**
   * \returns the name of the subclassed application
   */ 
  virtual std::string getName()=0;

  /**
   * \returns the description of the subclassed application
   */
  virtual std::string getDescription()=0;


  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

  /**
   * \brief Perform the set-up steps common to all crux commands:
   * initialize parameters, parse command line, set verbosity, open
   * output directory, write params file. 
   *
   * Uses the given command name, arguments and options for parsing the
   * command line.
   */
  virtual void initialize(
    const char** argument_list, ///< list of required arguments
    int num_arguments,          ///< number of elements in arguments_list
    const char** option_list,   ///< list of optional flags
    int num_options,            ///< number of elements in options_list
    int argc,                   ///< number of tokens on cmd line
    char** argv                 ///< array of command line tokens
  );


  /**
   * Frees an allocated CruxApplication
   */
  virtual ~CruxApplication();

  /**
   * Should this application be kept from the usage statement?
   */
  virtual bool hidden();

  /**
   * Writes the parameter file
   */
  virtual void writeParamFile();
};



#endif
