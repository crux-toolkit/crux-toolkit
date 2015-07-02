/**
 * \file CruxHardklorApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef CRUXHARDKLORAPPLICATION_H
#define CRUXHARDKLORAPPLICATION_H

#include "app/CruxApplication.h"

#include <string>
#include <fstream>

class CruxHardklorApplication: public CruxApplication {

 public:
  CruxHardklorApplication();
  ~CruxHardklorApplication();

  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CruxHardklorApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for CruxHardklorApplication
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

  /**
   * \brief runs hardklor on the input spectra
   * \returns whether hardklor was successful or not
   */
  static int main(
    const std::string& ms1 ///< file path of spectra to process
  );
  
 protected:
  static void addArg(
    std::vector<char*>* args,
    const std::string& arg
  );

  static void addArg(
    std::vector<char*>* args,
    const std::string& name,
    const std::string& value
  );

  static void addArg(
    std::vector<char*>* args,
    const std::string& name,
    bool value
  );
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
