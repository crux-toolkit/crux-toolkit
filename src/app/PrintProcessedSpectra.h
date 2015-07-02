/**
 * \file PrintProcessedSpectra.h
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 * REVISION:
 */
#ifndef PRINTPROCESSEDSPECTRA_H
#define PRINTPROCESSEDSPECTRA_H
#include "CruxApplication.h"


#include "util/crux-utils.h"
#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumCollection.h"
#include "model/FilteredSpectrumChargeIterator.h"
#include "model/Scorer.h"

#include <string>

class PrintProcessedSpectra: public CruxApplication {

 public:
  /**
   * \returns a blank PrintProcessedSpectra object
   */
  PrintProcessedSpectra();
  
  /**
   * Destructor
   */
  ~PrintProcessedSpectra();

  /**
   * main method for PrintProcessedSpectra
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns the command name for PrintProcessedSpectra
   */
  virtual std::string getName() const;

  /**
   * \returns the description for PrintProcessedSpectra
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
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem() const;

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


