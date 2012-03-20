/**
 * \file CruxBullseyeApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 4 November 2011
 * \brief Interface for calling hardklor.
 *****************************************************************************/
#ifndef CRUXBULLSEYEAPPLICATION_H
#define CRUXBULLSEYEAPPLICATION_H

#include "CruxApplication.h"

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

  /**
   * \returns the command name for CruxBullseyeApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for CruxBullseyeApplication
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
