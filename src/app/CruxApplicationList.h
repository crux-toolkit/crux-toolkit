/**
 * \file CruxApplicationList.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Maintains a list of executable applications.  Calling the main method
 * will try to find the application by name and try to execute it if found.
 *****************************************************************************/
#ifndef CRUXAPPLICATIONLIST_H
#define CRUXAPPLICATIONLIST_H

#include "CruxApplication.h"

#include <map>
#include <vector>
#include <string>


class CruxApplicationList {

 protected:
  std::vector<CruxApplication*> applications_; ///< list of applications
  std::map<int, std::string> messages_; ///< messages and their indexes
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
  void add(
    CruxApplication* application ///< application to add
  );

  /**
   * Adds a message to be printed in the usage statement
   */
  void addMessage(
    const std::string& message ///< message to add
  );

  /**
   * \returns an application by a name,
   * returns NULL if not found
   */
  CruxApplication* find(
    const std::string& appname ///< name of application to find
  );

  /**
   * \returns an application by a name,
   * returns NULL if not found
   */
  CruxApplication* find(
    const char* appname ///< name of application to find
  );

  /**
   * prints out the usage statement for this application list.
   * Each applications name is printed along with its description
   */
  void usage();

  /**
   * Gets the name of the list
   */
  std::string getListName() const;

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
