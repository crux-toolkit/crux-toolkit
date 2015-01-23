/**
 * \file PercolatorApplication.h 
 * AUTHOR: Manije Naseri
 * CREATE DATE: 4 September 2012
 * \
 *****************************************************************************/
#ifndef MAKEPINAPPLICATION_H
#define MAKEPINAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>
#include <string>

using namespace std;

class MakePinApplication: public CruxApplication {

 public:

  /**
   * \returns a blank MakePinApplication object
   */
  MakePinApplication();

  /**
   * Destructor
   */
  ~MakePinApplication();

  /**
   * main method for MakePinApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * runs make-pin application
   */
  static int main(std::vector<std::string>& paths);

  /**
   * \returns the command name for MakePinApplication
   */
  virtual std::string getName();

  /**
   * \returns the description for MakePinApplication
   */
  virtual std::string getDescription();

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory();

  
  /**
   * hide sequest search 
  */

  virtual bool hidden();
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
