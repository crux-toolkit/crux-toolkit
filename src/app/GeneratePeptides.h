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
#include "util/crux-utils.h"
#include "io/carp.h"
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
  virtual std::string getName() const;

  /**
   * \returns The description for GeneratePeptides.
   */
  virtual std::string getDescription() const;

  /**
   * \returns The command arguments
   */
  virtual std::vector<std::string> getArgs() const;

  /**
   * \returns The command options
   */
  virtual std::vector<std::string> getOptions() const;

  /**
   * \returns The command outputs
   */
  virtual std::vector< std::pair<std::string, std::string> > getOutputs() const;

  /**
   * \returns The file stem of the application, default getName.
   */
  virtual std::string getFileStem() const;

  /**
   * \returns The enum of the application, default MISC_COMMAND.
   */
  virtual COMMAND_T getCommand() const;

  /**
   * \returns False, i.e. GeneratePeptides does not require an
   * output directory.
   */
  virtual bool needsOutputDirectory() const;

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


