#ifndef TIDEINDEXAPPLICATION_H
#define TIDEINDEXAPPLICATION_H

#include "CruxApplication.h"

using namespace std;

class TideIndexApplication : public CruxApplication {

protected:

public:

  /**
   * Constructor
   */
  TideIndexApplication();

  /**
   * Destructor
   */
  ~TideIndexApplication();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  /**
   * Returns the command name
   */
  virtual string getName();

  /**
   * Returns the command description
   */
  virtual string getDescription();

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
