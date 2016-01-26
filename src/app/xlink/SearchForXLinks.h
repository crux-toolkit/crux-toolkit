/**
 * \file SearchForXLinks.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 6 December 2010
 * \brief Object for running search-for-xlinks
 *****************************************************************************/

#ifndef SEARCHFORXLINKS_H
#define SEARCHFORXLINKS_H

#include "app/CruxApplication.h"
#include "io/DelimitedFileReader.h"

#include <string>

class SearchForXLinks: public CruxApplication {

 protected:
  /**
   * The main method for the new xlink search code
   * \returns return code after execution
   */
  int xlinkSearchMain();

  /**
   * The main method for the old xlink search code
   * \returns return code after execution
   */
  int xhhcSearchMain();

 public:

  /**
   * \returns a blank SearchForXLinks object
   */
  SearchForXLinks();
  
  /**
   * Destructor
   */
  virtual ~SearchForXLinks();

  /**
   * main method for SearchForXLinks
   */
  virtual int main(int argc, char** argv);

  virtual void processParams();

  /**
   * \returns the command name for SearchForXLinks
   */
  virtual std::string getName() const;

  /**
   * \returns the description for SearchForXLinks
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
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand() const;

  virtual bool needsOutputDirectory() const;


};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
