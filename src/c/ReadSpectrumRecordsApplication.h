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

protected:

  void show(HeadedRecordReader& reader);

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

  virtual bool hidden();
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
