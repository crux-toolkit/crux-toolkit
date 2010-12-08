/**
 * \file SequestSearch.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running sequest-search
 *****************************************************************************/
#ifndef SEQUESTSEARCH_H
#define SEQUESTSEARCH_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"
#include "OutputFiles.h"


#include <string>
#include <vector>

class SequestSearch: public CruxApplication {

 protected:
 public:

  SequestSearch();
  ~SequestSearch();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();
  virtual std::string getFileStem();
};


#endif
