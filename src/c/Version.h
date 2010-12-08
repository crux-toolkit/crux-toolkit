/**
 * \file Version.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#ifndef VERSION_H
#define VERSION_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class Version: public CruxApplication {

 public:

  Version();
  ~Version();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
