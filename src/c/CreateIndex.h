/**
 * \file CreateIndex.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running create-index
 *****************************************************************************/
#ifndef CREATEINDEX_H
#define CREATEINDEX_H

#include "Index.h"
#include "parameter.h"
#include "CruxApplication.h"

#include <string>

class CreateIndex: public CruxApplication {

 public:
  /**
   * \returns a blank CreateIndex object
   */
  CreateIndex();
  
  /**
   * Destructor
   */
  ~CreateIndex();

  /**
   * main method for CreateIndex
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CreateIndex
   */
  virtual std::string getName();

  /**
   * \returns the description for CreateIndex
   */
  virtual std::string getDescription();

  /**
   * \returns the file stem of the application, default getName.
  */
  virtual std::string getFileStem();


  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  virtual bool hidden();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
