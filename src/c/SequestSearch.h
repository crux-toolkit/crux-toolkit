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

 public:
  /**
   * \returns a blank SequestSearch object
   */
  SequestSearch();
  
  /**
   * Destructor
   */
  ~SequestSearch();

  /**
   * main method for SequestSearch
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for SequestSearch
   */
  virtual std::string getName();

  /**
   * \returns the description for SequestSearch
   */
  virtual std::string getDescription();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
