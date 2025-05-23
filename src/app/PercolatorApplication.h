/**
 * \file PercolatorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef PERCOLATORAPPLICATION_H
#define PERCOLATORAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>


class PercolatorApplication: public CruxApplication {

 protected:

  //Calls the main method in Percolator Application
  static int percolatorMain(int argc, char* argv[]);


 public:

  /**
   * \returns a blank PercolatorApplication object
   */
  PercolatorApplication();

  /**
   * Destructor
   */
  ~PercolatorApplication();

  /**
   * main method for PercolatorApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for PercolatorApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for PercolatorApplication
   */
  virtual std::string getDescription() const;

  /**
   * \returns the command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns the command options
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * \returns the command outputs
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

  /**
   * \brief runs hardklor on the input spectra
   * \returns whether hardklor was successful or not
   */
  int main(
    const std::string& input_pinxml, // file path of spectra to process
    const std::string& output_dir_to_overwrite = "", // added by Yang
    const std::string& app_name = "" // added by Andy
  );
  
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
