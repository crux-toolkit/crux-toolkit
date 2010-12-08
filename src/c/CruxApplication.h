/**
 * \file CruxApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#ifndef CRUXAPPLICATION_H
#define CRUXAPPLICATION_H

#include <string>

class CruxApplication{
 public:
  /**
   * the main method for CruxApplication.  Subclasses of
   * CruxApplication define this.
   * \returns exit code for the executed program.   
   */
  virtual int main(int argc, char** argv)=0;

  /**
   * \returns the name of the subclassed application
   */ 
  virtual std::string getName()=0;

  /**
   * \returns the file stem of the application, default blank.
   */
  virtual std::string getFileStem();

  /**
   * \returns the description of the subclassed application
   */
  virtual std::string getDescription()=0;

  /**
   * Frees an allocated CruxApplication
   */
  virtual ~CruxApplication();

};



#endif
