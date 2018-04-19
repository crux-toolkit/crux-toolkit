/**
 * \file PercolatorAdapter.h
 * $Revision$
 * \brief Converts Percolator result objects to Crux result objects.
 */

#ifndef PERCOLATORADAPTER_H_
#define PERCOLATORADAPTER_H_

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "Caller.h"
#include "model/MatchCollection.h"
#include "model/PeptideMatch.h"
#include "model/PostProcessProtein.h"
#include "model/ProteinMatchCollection.h"
#include "model/ProteinMatch.h"
#include "ProteinProbEstimator.h"
#include "model/SpectrumMatch.h"
#include "Scores.h"

/**
 * \brief: Converts Percolator results objects to Crux results objects.
 * Class inherits the Caller class from percolator, which runs the percolator
 * algorithm and prints out the results.
 * All of the collections can then be accessed via the ProteinMatchCollection object.
 * During the execution of Caller::run(), percolator discards the psms after
 * Calculating the peptide level statistics, hence the need to pull the data for the
 * psms before the removal occurs.
 */
class PercolatorAdapter : public Caller {

 public:

  /**
   * Constructor for PercolatorAdapter.
   */
  PercolatorAdapter();

  /**
   * Destructor for PercolatorAdapter
   */
  virtual ~PercolatorAdapter();

  void deleteCollections();

  /**
   * Adds PSM scores from Percolator objects into a ProteinMatchCollection
   */
  void processPsmScores(Scores& allScores);

  /**
   * Adds peptide scores from Percolator objects into a ProteinMatchCollection
   */
  void processPeptideScores(Scores& allScores);

  /**
   * Adds protein scores from Percolator objects into a ProteinMatchCollection
   */
  void processProteinScores(ProteinProbEstimator* protEstimator);

  /*
   *\returns the ProteinMatchCollection, to be called after Caller::run() is finished
   */
  ProteinMatchCollection* getProteinMatchCollection();

  /**
   *\returns the decoy ProteinMatchCollection, to be called after Caller::run() is finished
   */
  ProteinMatchCollection* getDecoyProteinMatchCollection();

  static int findFeatureIndex(std::string feature);
  static std::map<int, int> mapChargeFeatures(); // map index of feature -> charge

  static double unnormalize(
    const PSMDescription* psm,
    int featureIndex,
    double* normDiv = NULL,
    double* normSub = NULL
  );

  static void printScores(Scores* scores, int label, std::ostream& os);
  
 protected:
    
  ProteinMatchCollection* collection_; ///< Collection containing all of the psm, peptide, and protein results.
  ProteinMatchCollection* decoy_collection_;  ///< Decoy ProteinMatchCollection
  std::vector<MatchCollection*> match_collections_made_; ///< MatchCollections created
  std::vector<PostProcessProtein*> proteins_made_; ///< Proteins created
  
  /**
   * Given a Percolator psm_id in the form ".*_([0-9]+)_[^_]*",
   * parse the file_idx, scan#, charge, and rank
   * Returns true on success
   */
  static bool parsePSMId(
    const std::string& psm_id, ///<psmid to parse
    int& file_idx, ///< psm file idx
    int& charge ///< psm charge
  );

  /**
   * Compare two AbstractMatches by Percolator score
   */
  static bool comparePercolatorScores(
    AbstractMatch* lhs, ///< first match with Percolator score to compare
    AbstractMatch* rhs ///< second match with Percolator score to compare
  );

  /**
   * \returns a Crux peptide from the PSM
   */
  Crux::Peptide* extractPeptide(
    PSMDescription* psm, ///< psm
    int charge_state, ///< charge state
    bool is_decoy ///< is psm a decoy?
    );

  PostProcessProtein* makeProtein(const std::string& name);
};

#endif /* PERCOLATORADAPTER_H_ */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

