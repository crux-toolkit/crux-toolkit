/**
 * \file MonoLinkPeptide.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a mono-link peptide in an xlink search
 *****************************************************************************/
#ifndef MONOLINKPEPTIDE_H_
#define MONOLINKPEPTIDE_H_

#include "model/objects.h"
#include "util/utils.h"

#include "LinearPeptide.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class MonoLinkPeptide : public LinearPeptide {
 protected:

 public:
  /**
   * Default constructor
   */
  MonoLinkPeptide();
  
  /**
   * Constructor with sequence
   */
  MonoLinkPeptide(
    char* sequence_ ///< sequence string
  );
  
  /**
   * Constructor from a Crux peptide
   */
  MonoLinkPeptide(
    Crux::Peptide* peptide ///< peptide object
  );

  /**
   * Default destructor
   */
  virtual ~MonoLinkPeptide();
  
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
  

};

/**
 * Comparison function of two linear peptide masses
 */
bool compareMonoLinkPeptideMass(
  const MonoLinkPeptide& xpep1, ///< linear peptide 1
  const MonoLinkPeptide& xpep2  ///< linear peptide 2
  );

/**
 * Comparison functions for lower_bound and upper_bound
 */
bool compareMonoLinkPeptideMassToFLOAT(
  const MonoLinkPeptide& xpep1, 
  FLOAT_T mass
  );


bool compareFLOATToMonoLinkPeptideMass(
  const FLOAT_T& mass,
  const MonoLinkPeptide& pep1);

#endif

/*
 * Local Variables:
 * mode: c 
 * c-basic-offset: 2
 * End:
 */
