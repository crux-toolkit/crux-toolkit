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
   * Perform Convert
   */
  virtual void convertFile(string input_format, string output_format, string input_file, string output_file_base, string database_file, bool distinct_matches);
  
  /**
   * Returns the command name
   */
  virtual string getName() const;

  /**
   * Returns the command description
   */
  virtual string getDescription() const;

  /**
   * Returns the command arguments
   */
  virtual vector<string> getArgs() const;

  /**
   * Returns the command options
   */
  virtual vector<string> getOptions() const;

  /**
   * \returns the outputs of the application as <name, description>
   */
  virtual vector< pair<string, string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;
  
};

#endif
