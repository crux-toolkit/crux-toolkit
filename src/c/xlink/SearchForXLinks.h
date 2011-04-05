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

  /**
   * \returns a blank SearchForXLinks object
   */
  SearchForXLinks();
  
  /**
   * Destructor
   */
  ~SearchForXLinks();

  /**
   * main method for SearchForXLinks
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for SearchForXLinks
   */
  virtual std::string getName();

  /**
   * \returns the description for SearchForXLinks
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
