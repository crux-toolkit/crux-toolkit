/**
 * \file LinearPeptide.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a Normal peptide in an xlink search
 *****************************************************************************/
#ifndef LINEARPEPTIDE_H_
#define LINEARPEPTIDE_H_

#include "model/objects.h"
#include "util/utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class LinearPeptide : public XLinkMatch {
 protected:

  Crux::Peptide* peptide_; ///< Peptide this linear peptide referes to
  char* sequence_; ///< sequence of the peptide
 public:
  /**
   * Default constructor
   */
  LinearPeptide();
  
  /**
   * Constructor with sequence
   */
  LinearPeptide(
    char* sequence_ ///< sequence string
  );
  
  /**
   * Constructor from a Crux peptide
   */
  LinearPeptide(
    Crux::Peptide* peptide ///< peptide object
  );

  /**
   * Default destructor
   */
  virtual ~LinearPeptide();
  
  /**
   *Add candidates to the XLinkMatchCollection that are linear
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass
    FLOAT_T max_mass, ///< max mass
    bool is_decoy,
    XLinkMatchCollection& candidates ///< Vector of possible candidate -inout
    );

  /**
   * returns the candidate type, either a deadlink or a linear candidate
   */
  virtual XLINKMATCH_TYPE_T getCandidateType();
  
  /**
   * \returns the sequence of the peptide in string format
   */
  virtual std::string getSequenceString();
 
  /**
   * \returns the mass of the peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
  );
  
  /**
   *\returns a shuffled version of the peptide
   */
  virtual void shuffle(
    std::vector<XLinkMatch*>& decoys
  );

  /**
   * predicts the ions for this peptide
   */
  virtual void predictIons(
    IonSeries* ion_series, ///< ion series to place the ions
    int charge ///< charge state of the peptide
    );
  
  /**
   * \returns the ion sequence as a string
   */
  std::string getIonSequence(
    Ion* ion ///< ion object
  );
  
  /**
   * \returns the peptide for this match
   */
  virtual Crux::Peptide* getPeptide( 
    int peptide_idx ///< should always be zero
  );
  
  /**
   * \returns the number of missed cleavages
   */
  virtual int getNumMissedCleavages();
  
  /**
   * \returns whether this peptide is modified by a variable mod
   */
  virtual bool isModified();

  /**
   *\returns "target" or "decoy"
   */
  string getDecoyType();

};

/**
 * Comparison function of two linear peptide masses
 */
bool compareLinearPeptideMass(
  const LinearPeptide& xpep1, ///< linear peptide 1
  const LinearPeptide& xpep2  ///< linear peptide 2
  );

/**
 * Comparison functions for lower_bound and upper_bound
 */
bool compareLinearPeptideMassToFLOAT(
  const LinearPeptide& xpep1, 
  FLOAT_T mass
  );


bool compareFLOATToLinearPeptideMass(
  const FLOAT_T& mass,
  const LinearPeptide& pep1);

#endif

/*
 * Local Variables:
 * mode: c 
 * c-basic-offset: 2
 * End:
 */
