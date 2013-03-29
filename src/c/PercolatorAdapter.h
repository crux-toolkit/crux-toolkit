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

#include "MatchCollection.h"
#include "PeptideMatch.h"
#include "PostProcessProtein.h"
#include "ProteinMatchCollection.h"
#include "ProteinMatch.h"
#include "SpectrumMatch.h"
#include "Scores.h"

using namespace std;

class PercolatorAdapter {

public:

  /**
   * Constructor for PercolatorAdapter. This should not be called, since all of
   * this class's functions are static.
   */
  PercolatorAdapter();

  /**
   * Destructor for PercolatorAdapter
   */
  virtual ~PercolatorAdapter();

  /**
   * Converts a set of Percolator scores into a Crux MatchCollection
   */
  static MatchCollection* psmScoresToMatchCollection(
    Scores* scores ///< percolator scores to convert
  );

  /**
   * Adds PSM scores from Percolator objects into a ProteinMatchCollection
   */
  static void addPsmScores(
    ProteinMatchCollection* collection, ///< collection to add scores to
    Scores* scores ///< percolator scores to add
  );

  /**
   * Adds protein scores from Percolator objects into a ProteinMatchCollection
   */
  static void addProteinScores(
    ProteinMatchCollection* collection, ///< collection to add scores to
    Scores* scores ///< percolator scores to add
  );

  /**
   * Adds peptide scores from Percolator objects into a ProteinMatchCollection
   */
  static void addPeptideScores(
    ProteinMatchCollection* collection, ///< collection to add scores to
    Scores* scores ///< percolator scores to add
  );

protected:
  /**
   * Given a Percolator psm_id in the form ".*_([0-9]+)_[^_]*",
   * find the charge state (matching group)
   */
  static int parseChargeState(
    string psm_id ///< psm to parse charge state from
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
  static Crux::Peptide* extractPeptide(
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
    string& seq, ///< sequence -out
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

