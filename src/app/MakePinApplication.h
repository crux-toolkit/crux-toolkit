/**
 * \file PercolatorApplication.h 
 * AUTHOR: Manije Naseri
 * CREATE DATE: 4 September 2012
 * \
 *****************************************************************************/
#ifndef MAKEPINAPPLICATION_H
#define MAKEPINAPPLICATION_H

#include "CruxApplication.h"

#include <string>
#include <fstream>
#include <string>

using namespace std;

class MakePinApplication: public CruxApplication {

 public:

  /**
   * \returns a blank MakePinApplication object
   */
  MakePinApplication();

  /**
   * Destructor
   */
  ~MakePinApplication();

  /**
   * main method for MakePinApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * runs make-pin application
   */
  static int main(const std::vector<std::string>& paths);

  /**
   * \returns the command name for MakePinApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for MakePinApplication
   */
  virtual std::string getDescription() const;

  /**
   * \returns the command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns the command options
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * \returns the command outputs
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory() const;

  
  /**
   * hide sequest search 
  */

  virtual bool hidden() const;
};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
