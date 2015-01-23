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

  /**
   * \returns the command name for CometApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for CometApplication
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

  /**
   * Sets the parameters for the Comet application using the crux parameters
   */
  void setCometParameters(
    std::vector<InputFileInfo*> &pvInputFiles, ///<vector of input spectra files
    CometSearchManager& searchMgr ///< SearchManager to set the parameters
    );
  
  /**
   * Extracts the variable modification info from the string
   */
  void calcVarMods(
    const char* var_str, ///< variable modification string
    VarMods& varmods ///< Variable modification structure to set
    );

  /**
   * Sets a double range from a string
   */
  void getDoubleRange(
    const char* str, ///< range string to set -in
    DoubleRange& doubleRangeParam ///< DoubleRange parameter -out
    );

  /*
   * sets the integer range from a string
   */
  void getIntRange(
    const char* str, ///< range string to set -in
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
