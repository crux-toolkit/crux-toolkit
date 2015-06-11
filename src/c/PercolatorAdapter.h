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
#include "MatchCollection.h"
#include "PeptideMatch.h"
#include "PostProcessProtein.h"
#include "ProteinMatchCollection.h"
#include "ProteinMatch.h"
#include "ProteinProbEstimator.h"
#include "SpectrumMatch.h"
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

  /**
   * Converts a set of Percolator scores into a Crux MatchCollection
   */
  void psmScoresToMatchCollection(
    MatchCollection** match_collection,  ///< out parameter for targets
    MatchCollection** decoy_match_collection ///< out parameter for decoys
  );

  /**
   * Adds PSM scores from Percolator objects into a ProteinMatchCollection
   */
  void addPsmScores();

  /**
   * Adds protein scores from Percolator objects into a ProteinMatchCollection
   */
  void addProteinScores();

  /**
   * Adds peptide scores from Percolator objects into a ProteinMatchCollection
   */
  void addPeptideScores();

  /*
   *\returns the ProteinMatchCollection, to be called after Caller::run() is finished
   */
  ProteinMatchCollection* getProteinMatchCollection();

  /**
   *\returns the decoy ProteinMatchCollection, to be called after Caller::run() is finished
   */
  ProteinMatchCollection* getDecoyProteinMatchCollection();

  int run();
  
protected:
    
  ProteinMatchCollection* collection_; ///< Collection containing all of the psm, peptide, and protein results.
  ProteinMatchCollection* decoy_collection_;  ///< Decoy ProteinMatchCollection
  std::vector<MatchCollection*> match_collections_made_; ///< MatchCollections created
  std::vector<PostProcessProtein*> proteins_made_; ///< Proteins created
  
  /**
   * Given a Percolator psm_id in the form ".*_([0-9]+)_[^_]*",
   * parse the file_idx, scan#, charge, and rank
   */
  static void parsePSMId(
    const std::string& psm_id, ///<psmid to parse
    int& file_idx, ///< psm file idx
    int& scan, ///< psm scan
    int& charge, ///< psm charge
    int& rank ///< psm rank
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

  /**
   * \returns the modified and unmodified peptide sequence
   * for the psm
   */
  static MODIFIED_AA_T* getModifiedAASequence(
    PSMDescription* psm, ///< psm -in
    std::string& seq, ///< sequence -out
    FLOAT_T& peptide_mass ///< calculated mass of peptide with modifications -out
    );
};

#endif /* PERCOLATORADAPTER_H_ */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

