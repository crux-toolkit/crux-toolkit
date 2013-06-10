#ifndef TIDERESULTSAPPLICATION_H
#define TIDERESULTSAPPLICATION_H

#include "CruxApplication.h"

#include <fstream>
#include "tide/abspath.h"
#include "tide/results_writer.h"

using namespace std;

class TideResultsApplication : public CruxApplication {

protected:

public:

  /**
   * Constructor
   */
  TideResultsApplication();

  /**
   * Destructor
   */
  ~TideResultsApplication();

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
