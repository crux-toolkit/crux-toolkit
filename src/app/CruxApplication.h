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
#include <vector>
#include <map>

class CruxApplication{
 public:
  /**
   * the main method for CruxApplication.  Subclasses of
   * CruxApplication define this.
   * \returns exit code for the executed program.   
   */
  virtual int main(int argc, char** argv) = 0;

  /**
   * \returns the name of the subclassed application
   */ 
  virtual std::string getName() const = 0;

  /**
   * \returns the description of the subclassed application
   */
  virtual std::string getDescription() const = 0;

  /**
   * \returns the arguments of the application
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns the options of the application
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * \returns the outputs of the application as name -> description
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem() const;

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand() const;

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory() const;

  virtual void initialize(int argc, char** argv);

  /**
   * Frees an allocated CruxApplication
   */
  virtual ~CruxApplication();

  /**
   * Should this application be kept from the usage statement?
   */
  virtual bool hidden() const;

  /**
   * Read in all parameters from command line and parameter file
   */
  static void initializeParams(
    const std::string& appName,
    const std::vector<std::string>& appArgs,
    const std::vector<std::string>& appOptions,
    int argc,
    char** argv
  );

  /**
   * Process parameters after they have been set up, but before they have been
   * finalized
   */
  virtual void processParams();

  /**
   * Get usage statement for the program.
   */
  static std::string getUsage(
    const std::string& appName,
    const std::vector<std::string>& args,
    const std::vector<std::string>& options
  );
};



#endif
