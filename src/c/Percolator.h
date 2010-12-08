/**
 * \file Percolator.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running percolator
 *****************************************************************************/
#ifndef PERCOLATOR_H
#define PERCOLATOR_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class Percolator: public CruxApplication {

 public:

  Percolator();
  ~Percolator();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
