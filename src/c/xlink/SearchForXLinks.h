/**
 * \file SearchForXLinks.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef SEARCHFORXLINKS_H
#define SEARCHFORXLINKS_H

#include "CruxApplication.h"
#include "DelimitedFileReader.h"

#include <string>

class SearchForXLinks: public CruxApplication {

 public:

  SearchForXLinks();
  ~SearchForXLinks();
  virtual int main(int argc, char** argv);
  virtual std::string getName();
  virtual std::string getDescription();

};


#endif
