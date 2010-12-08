/**
 * \file CruxApplicationList.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for a CruxApplicationList
 *****************************************************************************/
#ifndef CRUXAPPLICATIONLIST_H
#define CRUXAPPLICATIONLIST_H

#include "CruxApplication.h"

#include <vector>
#include <string>

class CruxApplicationList : std::vector<CruxApplication*> {

 protected:
  std::string list_name_; ///<Name of this list

 public:
  /**
   * \brief Creates an application list with a listname
   */
  CruxApplicationList(const char* list_name);

  
  /**
   * Frees an allocated CruxApplicationList
   */
  ~CruxApplicationList();

  /**
   * Adds an application pointer to the list of applications
   */
  void add(CruxApplication* application);

  /**
   * \returns an application by a name,
   * returns NULL if not found
   */
  CruxApplication* find(const std::string& appname);

  /**
   * \returns an application by a name,
   * returns NULL if not found
   */
  CruxApplication* find(const char* appname);

  /**
   * prints out the usage statement for this application list.
   * Each applications name is printed along with its description
   */
  void usage();

  /**
   * the main method for CruxApplicationList.  Attempts to find
   * an application by name from the first argument.  If successful,
   * calls that applications main method with the rest of the parameters.
   */
  int main(int argc, char** argv);
};


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
