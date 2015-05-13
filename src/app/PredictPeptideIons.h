/**
 * \file PredictPeptideIons.h
 *
 * AUTHOR: Manijeh Naseri
 * CREATE DATE: January 27, 2012
 * DESCRIPTION: Main method for the predict-peptide-ions.
 *              Given a peptide sequence, and a charge state, predict
 *              the fragmentation ions.
 */
#ifndef PREDICTPEPTIDEIONS_H
#define PREDICTPEPTIDEIONS_H

#include "CruxApplication.h"
#include "util/crux-utils.h"
#include "io/carp.h"
#include "parameter.h"

class PredictPeptideIons: public CruxApplication {

 public:
  /**
   * \returns A blank PredictPeptideIons object.
   */
  PredictPeptideIons();
  
  /**
   * Destructor
   */
  ~PredictPeptideIons();

  /**
   * Main method for PredictPeptideIons.
   */
  virtual int main(int argc, char** argv);

  /**
   * \returns The command name for PredictPeptideIons.
   */
  virtual std::string getName() const;

  /**
   * \returns The description for PredictPeptideIons.
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
   * \returns The enum of the application, PREDICT_PEPTIDE_IONS_COMMAND.
   */
  virtual COMMAND_T getCommand() const;

};


#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


