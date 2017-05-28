#ifndef READSPECTRUMRECORDSAPPLICATION_H
#define READSPECTRUMRECORDSAPPLICATION_H

#include "CruxApplication.h"

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <gflags/gflags.h>
#include "tide/spectrum_collection.h"
#include "tide/records.h"

using namespace std;

class ReadSpectrumRecordsApplication : public CruxApplication {

 public:

  /**
   * Constructor
   */
  ReadSpectrumRecordsApplication();

  /**
   * Destructor
   */
  ~ReadSpectrumRecordsApplication();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  /**
   * Returns the command name
   */
  virtual string getName() const;

  /**
   * Returns the command description
   */
  virtual string getDescription() const;

  virtual vector<string> getArgs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

  virtual bool hidden() const;
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
