/**
 * \file PrintVersion.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#ifndef PRINTVERSION_H
#define PRINTVERSION_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class PrintVersion: public CruxApplication {

 public:

  PrintVersion();
  ~PrintVersion();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
