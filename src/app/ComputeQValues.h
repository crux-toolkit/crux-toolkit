/**
 * \file ComputeQValues.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running compute-q-values
 *****************************************************************************/
#ifndef ComputeQValues_H
#define ComputeQValues_H

#include "CruxApplication.h"

#include <string>

class ComputeQValues: public CruxApplication {

 public:

  /**
   * \returns a blank ComputeQValues object
   */
  ComputeQValues();
  
  /**
   * Destructor
   */
  ~ComputeQValues();

  /**
   * main method for ComputeQValues
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for ComputeQValues
   */
  virtual std::string getName();

  /**
   * \returns the description for ComputeQValues
   */
  virtual std::string getDescription();

  /**
   * \returns the filestem for ComputeQValues
   */
  virtual std::string getFileStem();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns whether the application needs the output directory or not.
   */
  virtual bool needsOutputDirectory();


};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
