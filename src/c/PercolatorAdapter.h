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
#include "Caller.h"
using namespace std;

/**
 * \brief: Converts Percolator results objects to Crux results objects.
 * Class inherits the Caller class from percolator, which runs the percolator
 * algorithm and prints out the results.  During the printing of the results, the
 * PercolatorAdapter class overrides the methods for printing out the xml result file
 * to collect the psm, peptide, and protein results from the internal fullset.
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
  
protected:
    
  ProteinMatchCollection* collection_; ///< Collection containing all of the psm, peptide, and protein results.
  ProteinMatchCollection* decoy_collection_;  ///< Decoy ProteinMatchCollection
  vector<MatchCollection*> match_collections_made_; ///< MatchCollections created
  vector<PostProcessProtein*> proteins_made_; ///< Proteins created
  
  /**
   * Calls Percolator's overridden Caller::writeXML_PSMs() and then
   * Collects all of the psm results
   */
  virtual void writeXML_PSMs();
  
  /**
   * Calls Percolator's overridden Caller::writeXML_Peptides() and then
   * Collects all of the peptide results from the fullset
   */
  virtual void writeXML_Peptides();
  
  /**
   * Calls Percolator's overriden Caller::writeXMLProteins() and then
   * Collects all of the protein results from fido
   */
  virtual void writeXML_Proteins();
  
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

