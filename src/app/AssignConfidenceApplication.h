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
#include "model/objects.h"
#include "parameter.h"
#include "model/Protein.h"
#include "model/Spectrum.h"
#include "io/SpectrumCollection.h"
#include "model/Scorer.h"
#include "model/Match.h"
#include "model/MatchCollection.h"
#include "io/OutputFiles.h"
#include "model/Peptide.h"
#include "boost/tuple/tuple.hpp" // This will be <tuple> once we move to C++11.
#include "boost/tuple/tuple_comparison.hpp"

/**
 * Legal values for the --estimation-method option.
 */
enum _estimation_method {
  INVALID_METHOD,
  MIXMAX_METHOD,       ///< compute q-values using mix-max (Uri Keich)
  TDC_METHOD,          ///< compute q-values using TDC (Elias-Gygi)  
  PEPTIDE_LEVEL_METHOD,///< simple peptide-level decoy-based estimation (WOTE)
  NUMBER_METHOD_TYPES  ///< always keep this last so the value
                        /// changes as cmds are added
};

typedef enum _estimation_method ESTIMATION_METHOD_T;

class AssignConfidenceApplication : public CruxApplication {
 protected:
  map<pair<string, unsigned int>, bool>* spectrum_flag_;  // this variable is used in Cascade Search, this is an idicator 
  unsigned int iteration_cnt_;
  OutputFiles* output_;
  unsigned int accepted_psms_;
  string index_name_;
  bool is_final_;

  class AtdcScoreSet {
   public:
    AtdcScoreSet(
      const std::vector<FLOAT_T>& targetScores,
      const std::vector< std::vector<FLOAT_T> >& decoyScores,
      bool ascending);
    static void getScores(
      MatchCollection* targets,
      std::map<int, MatchCollection*> decoys,
      SCORER_TYPE_T scoreType,
      bool ascending,
      std::vector<FLOAT_T>& outTargetScores,
      std::vector< std::vector<FLOAT_T> >& outDecoyScores);
    std::vector<FLOAT_T> fdps() const;
   private:
    void histBin(std::vector<int>& hist, FLOAT_T x) const;
    static bool sortScoresAsc(
      const boost::tuple<FLOAT_T, int, int, int>& x,
      const boost::tuple<FLOAT_T, int, int, int>& y);
    static bool sortScoresDesc(
      const boost::tuple<FLOAT_T, int, int, int>& x,
      const boost::tuple<FLOAT_T, int, int, int>& y);
    static bool sortScoresByIdentifier(
      const boost::tuple<FLOAT_T, int, int, int>& x,
      const boost::tuple<FLOAT_T, int, int, int>& y);

    std::vector< std::pair<FLOAT_T, std::vector<FLOAT_T> > > scores_; // <target score, [decoy scores]>
  };

 public:
  map<pair<string, unsigned int>, bool>* getSpectrumFlag();
  void setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag);
  void setIterationCnt(unsigned int iteration_cnt);
  void setOutput(OutputFiles *output);
  unsigned int getAcceptedPSMs();
  std::string getPeptideSeq(Crux::Match* match);

  /**
  * stores the name of the index file used in an iteration in Cascade Search.
  */
  void setIndexName(string index_name);

  /**
  * Stores a Boolean indicating whether this is the last iteration of the cascade.
  */
  void setFinalIteration(bool is_final);

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

  virtual int main(const vector<string>& input_files);

  static int getDirection(SCORER_TYPE_T scoreType);

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

  FLOAT_T* compute_qvalues_from_pvalues(
    FLOAT_T* pvalues,
    int      num_pvals,
    FLOAT_T  pi_zero);

  void peptide_level_filtering(
    MatchCollection* match_collection,
    std::map<string, FLOAT_T>* BestPeptideScore,
    SCORER_TYPE_T score_type,
    bool ascending);
  
  void identify_best_psm_per_peptide
    (MatchCollection* all_matches,
    SCORER_TYPE_T score_type);
  static void convert_fdr_to_qvalue(
    std::vector<FLOAT_T>& qvalues); ///< Come in as FDRs, go out as q-values.

  map<FLOAT_T, FLOAT_T> store_arrays_as_hash(
    const std::vector<FLOAT_T>& keys,
    const std::vector<FLOAT_T>& values);
  std::vector<FLOAT_T> compute_decoy_qvalues_tdc(
    std::vector<FLOAT_T>& target_scores,
    std::vector<FLOAT_T>& decoy_scores,
    bool ascending,
    FLOAT_T pi_zero);
  std::vector<FLOAT_T> compute_decoy_qvalues_mixmax(
    std::vector<FLOAT_T>& target_scores,
    std::vector<FLOAT_T>& decoy_scores,
    bool ascending,
    FLOAT_T pi_zero);
};

#endif //ASSIGNCONFIDENCE_H

