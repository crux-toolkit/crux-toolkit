#ifndef PSMCONVERTAPPLICATION_H
#define PSMCONVERTAPPLICATION_H

#include "CruxApplication.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include <string>

using namespace std;

class PSMConvertApplication : public CruxApplication {

public:

  /**
   * Constructor
   */
  PSMConvertApplication();

  /**
   * Destructor
   */
  ~PSMConvertApplication();

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

  /*
   * Accepts two strings, returns whether the first string ends with the second
   */
  virtual int endsWith(string s, string ending);

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();
  
};

#endif
