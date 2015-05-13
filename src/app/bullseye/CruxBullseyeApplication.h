/**
 * \file CruxBullseyeApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef CRUXBULLSEYEAPPLICATION_H
#define CRUXBULLSEYEAPPLICATION_H

#include "app/CruxApplication.h"

#include <string>
#include <fstream>

class CruxBullseyeApplication: public CruxApplication {

 protected:

  //Calls the main method in bullseye
  int bullseyeMain(int argc, char* argv[]);

 public:

  /**
   * \returns a blank CruxBullseyeApplication object
   */
  CruxBullseyeApplication();

  /**
   * Destructor
   */
  ~CruxBullseyeApplication();

  /**
   * main method for CruxBullseyeApplication
   */
  virtual int main(int argc, char** argv);

  int main(const std::string& input_ms1,
           const std::string& input_ms2,
           const std::string& match_ms2,
           const std::string& nomatch_ms2);

  /**
   * \returns the command name for CruxBullseyeApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for CruxBullseyeApplication
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
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
