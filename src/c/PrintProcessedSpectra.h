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


#include "crux-utils.h"
#include "carp.h"
#include "parameter.h"
#include "SpectrumCollection.h"
#include "FilteredSpectrumChargeIterator.h"
#include "Scorer.h"

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
  virtual std::string getName();

  /**
   * \returns the description for PrintProcessedSpectra
   */
  virtual std::string getDescription();

  /**
   * \returns the file stem of the application, default getName.
   */
  virtual std::string getFileStem();

  /**
   * \returns the enum of the application, default MISC_COMMAND
   */
  virtual COMMAND_T getCommand();

  virtual bool needsOutputDirectory();


};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


