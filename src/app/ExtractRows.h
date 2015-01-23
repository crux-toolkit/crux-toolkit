/**
 * \file ExtractRows.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Given a tab delimited file and column name and value, print
 * out all rows that pass the relation operator (default equals).
 *****************************************************************************/
#ifndef EXTRACTROWS_H
#define EXTRACTROWS_H

#include "CruxApplication.h"

#include <string>

class ExtractRows: public CruxApplication {

 public:

  /**
   * \returns a blank ExtractRows object
   */
  ExtractRows();

  /**
   * Destructor
   */
  ~ExtractRows();

  /**
   * main method for ExtractRows
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for ExtractRows
   */
  virtual std::string getName();

  /**
   * \returns the description for ExtractRows
   */
  virtual std::string getDescription();

};


#endif //EXTRACTROWS_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

