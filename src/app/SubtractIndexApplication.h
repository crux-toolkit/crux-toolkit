/**
 * \file CascadeSearchApplication.h 
 * AUTHOR: Attila Kertesz-Farkas
 * CREATE DATE: 30 May 2015
 * \brief Iterative PSM meta-search via Cascade protocol
 ***********************************************************/
#ifndef SUBTRACTINDEXAPPLICATION_H
#define SUBTRACTINDEXAPPLICATION_H

#include "CruxApplication.h"
#include "peptides.pb.h"
#include "TideSearchApplication.h"
#include <string>



class SubtractIndexApplication: public CruxApplication {

 public:

  /**
   * \returns a blank CascadeSearchApplication object
   */
  SubtractIndexApplication();
  
  /**
   * Destructor
   */
  ~SubtractIndexApplication();

  /**
   * main method for CascadeSearchApplication
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for CascadeSearchApplication
   */
  virtual std::string getName() const;

  /**
   * \returns the description for CascadeSearchApplication
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
   * \returns the filestem for CascadeSearchApplication
   */
  virtual std::string getFileStem() const;

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand() const;

  /**
   * \returns whether the application needs the output directory or not.
   */
  virtual bool needsOutputDirectory() const;

  virtual void processParams();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
