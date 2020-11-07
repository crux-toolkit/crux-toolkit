/**
 * \file KojakApplication.h 
 * AUTHOR: Michael Hoopmann
 * CREATE DATE: 11 October 2018
 * \brief Interface for calling Kojak.
 *****************************************************************************/
#ifndef KOJAKAPPLICATION_H
#define KOJAKAPPLICATION_H

#include "CruxApplication.h"
#include "KojakManager.h"
#include <string>
#include <fstream>
#include <iostream>

class KojakApplication: public CruxApplication {
public:

  /**
  * \Default constructor
  */
  KojakApplication();

  /**
  * \Default destructor
  */
  ~KojakApplication();

  /**
  * \returns the command name for KojakApplication
  */
  virtual std::string getName() const;

  /**
  * \returns the description for KojakApplication
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

  /**
   * Run param-medic and sets the parameters accordingly.
   */
  virtual void processParams();

  /**
  * \brief runs kojak on the input spectra
  * \returns whether kojak was successful or not
  */
  virtual int main(int argc, char** argv);
  int main(const std::vector<std::string>& input_files);

protected:
  KojakManager searchManager_;
};


#endif

/*
* Local Variables:
* mode: c
* c-basic-offset: 2
* End:
*/
