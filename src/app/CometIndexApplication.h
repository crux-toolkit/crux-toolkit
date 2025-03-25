/**
 * \file CruxHardklorApplication.h
 * AUTHOR: Charles Grant
 * CREATE DATE: March 1, 2025
 * \brief Interface for calling comet fragment ion indexer.
 *****************************************************************************/
#ifndef COMETINDEXAPPLICATION_H
#define COMETINDEXAPPLICATION_H

#include "CometApplication.h"
#include "CometSearch/Common.h"
#include "CometSearch/CometSearchManager.h"
#include <string>
#include <fstream>
#include <iostream>

class CometIndexApplication : public CometApplication {

public:
  /**
   * \returns a blank CometApplication object
   */
  CometIndexApplication();

  /**
   * Destructor
   */
  ~CometIndexApplication();

  /**
   * main method for CometApplication
   */
  virtual int main(int argc, char **argv);

  int main(void);

  /**
   * Sets the parameters for the Comet application using the crux parameters
  */
  virtual void setCometParameters(
    const vector<string>& spec_files,
    vector<InputFileInfo*>& pvInputFiles ///<vector of input spectra files
  );
  /**
   * \returns the command name for CometIndexApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for CometIndexApplication
   */
  virtual std::string getDescription() const;

  /**
   * \returns the command options
   */
  vector<string> getOptions() const;

  /**
   * \returns the command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns the command outputs
   */
  virtual std::vector<std::pair<std::string, std::string>> getOutputs() const;

  virtual COMMAND_T getCommand() const;

  /**
   * Processes the parameters
   */
  virtual void processParams();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
