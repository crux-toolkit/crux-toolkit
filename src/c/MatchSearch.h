/**
 * \file MatchSearch.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-matches
 *****************************************************************************/
#ifndef MATCHSEARCH_H
#define MATCHSEARCH_H

#include "CruxApplication.h"

#include <string>

class MatchSearch : public CruxApplication {
 protected:

 public:
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getFileStem();
  virtual std::string getDescription();

  MatchSearch();
  virtual ~MatchSearch();

};



#endif
