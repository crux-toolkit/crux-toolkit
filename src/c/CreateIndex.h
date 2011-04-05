/**
 * \file CreateIndex.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running create-index
 *****************************************************************************/
#ifndef CREATEINDEX_H
#define CREATEINDEX_H

#include "index.h"
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
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
