#ifndef LIBRAAPPLICATION_H
#define LIBRAAPPLICATION_H

#include "CruxApplication.h"

#include <iostream>
#include <fstream>
#include <vector>


class LibraApplication : public CruxApplication {

 public:

  /**
   * Constructor
   */
  LibraApplication();

  /**
   * Destructor
   */
  ~LibraApplication();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  int main(const std::vector<std::string>& input_files);

  int main(const std::vector<std::string>& input_files, const std::string input_index);

  /**
   * Returns the command name
   */
  virtual std::string getName() const;

  /**
   * Returns the command description
   */
  virtual std::string getDescription() const;

  /**
   * Returns the command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * Returns the command options
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * Returns the command outputs
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
