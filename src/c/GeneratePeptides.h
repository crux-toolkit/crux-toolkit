/**
 * \file GeneratePeptides.h
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: January 20, 2012
 * DESCRIPTION: Main method for the generate-peptides command.
 *              Output all peptide sequences in the given fasta file
 *              that fall within all peptide constraints.
 */
#ifndef GENERATEPEPTIDES_H
#define GENERATEPEPTIDES_H

#include "CruxApplication.h"
#include "crux-utils.h"
#include "carp.h"
#include "parameter.h"

class GeneratePeptides: public CruxApplication {

 public:
  /**
   * \returns A blank GeneratePeptides object.
   */
  GeneratePeptides();
  
  /**
   * Destructor
   */
  ~GeneratePeptides();

  /**
   * Main method for GeneratePeptides.
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns The command name for GeneratePeptides.
   */
  virtual std::string getName();

  /**
   * \returns The description for GeneratePeptides.
   */
  virtual std::string getDescription();

  /**
   * \returns The file stem of the application, default getName.
   */
  virtual std::string getFileStem();

  /**
   * \returns The enum of the application, default MISC_COMMAND.
   */
  virtual COMMAND_T getCommand();

  /**
   * \returns False, i.e. GeneratePeptides does not require an
   * output directory.
   */
  virtual bool needsOutputDirectory();

 protected:
  /**
   * Print header lines to stdout.
   */
  void printHeader();

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


