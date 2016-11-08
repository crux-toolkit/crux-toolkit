/**
 * \file CruxApplication.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Abstract Object for a CruxApplication
 *****************************************************************************/
#ifndef CRUXAPPLICATION_H
#define CRUXAPPLICATION_H
#include "model/objects.h"

#include <algorithm>
#include <string>
#include <vector>
#include <map>

class CruxApplication{
 public:
  /**
   * the main method for CruxApplication.  Subclasses of
   * CruxApplication define this.
   * \returns exit code for the executed program.   
   */
  virtual int main(int argc, char** argv) = 0;

  /**
   * \returns the name of the subclassed application
   */ 
  virtual std::string getName() const = 0;

  /**
   * \returns the description of the subclassed application
   */
  virtual std::string getDescription() const = 0;

  /**
   * \returns the arguments of the application
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns the options of the application
   */
  virtual std::vector<std::string> getOptions() const;

  template<typename T>
  static void addOptionsFrom(std::vector<std::string>* optionsVector) {
    using std::find;
    using std::string;
    using std::vector;
    CruxApplication* app = new T();
    vector<string> options = app->getOptions();
    delete app;
    for (vector<string>::const_iterator i = options.begin(); i != options.end(); i++) {
      if (find(optionsVector->begin(), optionsVector->end(), *i) == optionsVector->end()) {
        optionsVector->push_back(*i);
      }
    }
  }

  static void removeOptionFrom(std::vector<std::string>* optionsVector, std::string option) {
    using std::find;
    using std::string;
    using std::vector;
    vector<string>::iterator i;
    i = find(optionsVector->begin(), optionsVector->end(), option);
    if (i != optionsVector->end()) {
      optionsVector->erase(i);
    }
  }

  /**
   * \returns the outputs of the application as <name, description>
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  template<typename T>
  static void addOutputsFrom(
    std::vector< std::pair<std::string, std::string> >* outputsVector) {
    using std::pair;
    using std::set;
    using std::string;
    using std::vector;
    set<string> existing;
    for (vector< pair<string, string> >::const_iterator i = outputsVector->begin();
         i != outputsVector->end();
         i++) {
      existing.insert(i->first);
    }

    CruxApplication* app = new T();
    vector< pair<string, string> > outputs = app->getOutputs();
    delete app;
    for (vector< pair<string, string> >::const_iterator i = outputs.begin();
         i != outputs.end();
         i++) {
      if (existing.insert(i->first).second) {
        outputsVector->push_back(*i);
      }
    }
  }

  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem() const;

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand() const;

  /**
   * \returns whether the application needs the output directory or not. (default false).
   */
  virtual bool needsOutputDirectory() const;

  virtual void initialize(int argc, char** argv);

  /**
   * Frees an allocated CruxApplication
   */
  virtual ~CruxApplication();

  /**
   * Should this application be kept from the usage statement?
   */
  virtual bool hidden() const;

  /**
   * Read in all parameters from command line and parameter file
   */
  static void initializeParams(
    const std::string& appName,
    const std::vector<std::string>& appArgs,
    const std::vector<std::string>& appOptions,
    int argc,
    char** argv
  );

  /**
   * Process parameters after they have been set up, but before they have been
   * finalized
   */
  virtual void processParams();

  /**
   * Get usage statement for the program.
   */
  static std::string getUsage(
    const std::string& appName,
    const std::vector<std::string>& args,
    const std::vector<std::string>& options,
    bool full
  );
};



#endif
