/**
* \file AssignConfidenceApplication.h
* AUTHOR: Attila Kertesz-Farkas
* CREATE DATE: May 01, 2015
****************************************************************************/

#ifndef ASSIGNCONFIDENCE_H
#define ASSIGNCONFIDENCE_H

#include <stdlib.h>
#include <stdio.h>
#include "io/carp.h"
#include "util/crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "model/Protein.h"
#include "model/Spectrum.h"
#include "io/SpectrumCollection.h"
#include "model/Scorer.h"
#include "model/Match.h"
#include "model/MatchCollection.h"
#include "io/OutputFiles.h"
#include "model/Peptide.h"

class AssignConfidenceApplication : public CruxApplication {
protected: 
  map<pair<string, unsigned int>, bool>* spectrum_flag_;  // this variable is used in Cascade Search, this is an idicator 
  double cascade_fdr_;
  unsigned int iteration_cnt_;
  OutputFiles* output_;
  unsigned int accepted_psms_;
  string index_name_;
  bool combine_modified_peptides_;
  bool combine_charge_states_;
public:

  map<pair<string, unsigned int>, bool>* getSpectrumFlag();
  void setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag);
  void setCascadeFDR(double cascade_fdr);
  void setIterationCnt(unsigned int iteration_cnt);
  void setOutput(OutputFiles *output);
  unsigned int getAcceptedPSMs();
  std::string getPeptideSeq(Crux::Match* match);


  /**
  * stores the name of the index file used in an iteration in Cascade Search.
  */
  void setIndexName(string index_name);

  /**
  * \returns a blank ComputeQValues object
  */
  AssignConfidenceApplication();

  /**
  * Destructor
  */
  ~AssignConfidenceApplication();

  /**
  * main method for ComputeQValues
  */
  virtual int main(int argc, char** argv);

  virtual int main(const vector<string> input_files);

  /**
  * \returns the command name for ComputeQValues
  */
  virtual std::string getName() const;

  /**
  * \returns the description for ComputeQValues
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
  virtual std::vector< std::pair<std::string, std::string> >  getOutputs() const;

  /**
  * \returns the filestem for ComputeQValues
  */
  virtual std::string getFileStem() const;

  /**
  * \returns the enum of the application, default MISC_COMMAND
  */
  virtual COMMAND_T getCommand() const;
  
  /**
  * \Preprocesses the program arguments before finalizing them.
  */
  virtual void processParams();

  /**
  * \returns whether the application needs the output directory or not.
  */
  virtual bool needsOutputDirectory() const;

  FLOAT_T* compute_decoy_qvalues_tdc(
    FLOAT_T* target_scores,
    int      num_targets,
    FLOAT_T* decoy_scores,
    int      num_decoys,
    bool     reverse);

  FLOAT_T* compute_qvalues_from_pvalues(
    FLOAT_T* pvalues,
    int      num_pvals,
    FLOAT_T  pi_zero);

  FLOAT_T estimate_pi0(FLOAT_T* target_scores,
    int      num_targets,
    FLOAT_T* decoy_scores,
    int      num_decoys,
    bool     ascending);

  void peptide_level_filtering(
    MatchCollection* match_collection,
    std::map<string, FLOAT_T>* BestPeptideScore,
    SCORER_TYPE_T score_type,
    bool ascending);
  
  void identify_best_psm_per_peptide
    (MatchCollection* all_matches,
    SCORER_TYPE_T score_type);
  void convert_fdr_to_qvalue
    (FLOAT_T* qvalues,     ///< Come in as FDRs, go out as q-values.
    int      num_values);
  map<FLOAT_T, FLOAT_T>* store_arrays_as_hash
    (FLOAT_T* keys,
    FLOAT_T* values,
    int      num_values
    );
  FLOAT_T* compute_decoy_qvalues_tdc(
    FLOAT_T* target_scores,
    int      num_targets,
    FLOAT_T* decoy_scores,
    int      num_decoys,
    bool     forward,
    FLOAT_T  pi_zero
    );
  FLOAT_T* compute_decoy_qvalues_mixmax(
    FLOAT_T* target_scores,
    int      num_targets,
    FLOAT_T* decoy_scores,
    int      num_decoys,
    bool     ascending,
    FLOAT_T  pi_zero
    );

};


#endif //ASSIGNCONFIDENCE_H

