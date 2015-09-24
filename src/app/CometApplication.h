/**
 * \file CruxHardklorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef COMETAPPLICATION_H
#define COMETAPPLICATION_H

#include "CruxApplication.h"
#include "CometSearch/CometSearchManager.h"
#include <string>
#include <fstream>
#include <iostream>

class CometApplication: public CruxApplication {

 protected:

 public:

  /**
   * \returns a blank CometApplication object
   */
  CometApplication();

  /**
   * Destructor
   */
  ~CometApplication();

  /**
   * main method for CometApplication
   */
  virtual int main(int argc, char** argv);

  int main(const std::vector<std::string>& input_files);

  /**
   * \returns the command name for CometApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for CometApplication
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

  virtual COMMAND_T getCommand() const;

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory() const;

  virtual void processParams();

  /**
   * Sets the parameters for the Comet application using the crux parameters
   */
  void setCometParameters(
    const std::vector<std::string>& spec_files,
    std::vector<InputFileInfo*> &pvInputFiles, ///<vector of input spectra files
    CometSearchManager& searchMgr ///< SearchManager to set the parameters
    );
  
  /**
   * Extracts the variable modification info from the string
   */
  void calcVarMods(
    const std::string& var_str, ///< variable modification string
    VarMods& varmods ///< Variable modification structure to set
    );

  /**
   * Sets a double range from a string
   */
  void getDoubleRange(
    const std::string& str, ///< range string to set -in
    DoubleRange& doubleRangeParam ///< DoubleRange parameter -out
    );

  /*
   * sets the integer range from a string
   */
  void getIntRange(
    const std::string& str, ///< range string to set -in
    IntRange& intRangeParam ///< IntRange parameter -out
    );

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
