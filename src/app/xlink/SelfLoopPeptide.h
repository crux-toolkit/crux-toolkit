/**
 * \file SelfLoopPeptide.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a self-loop peptide in an xlink search
 *****************************************************************************/
#ifndef SELFLOOPPEPTIDE_H_
#define SELFLOOPPEPTIDE_H_

#include "model/objects.h"
#include "util/utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class SelfLoopPeptide : public XLinkMatch {
 protected:
  SelfLoopPeptide* target_;
  XLinkablePeptide linked_peptide_; ///< linkable peptide
  std::vector<int> link_pos_idx_; ///< the two link indices on the peptide
 public:
  
  /**
   * \returns the link position
   */
  int getLinkPos(
    int link_idx ///< link index (0 or 1)
  );

  /**
   * Default Constructor
   */
  SelfLoopPeptide();
  
  /**
   * Constructor that defines a linkable peptide and the positions
   */
  SelfLoopPeptide(
    XLinkablePeptide& peptide, ///< linkable peptide
    int posA, ///< 1st link position on peptide
    int posB  ///< 2nd link position on peptide
    );

  /**
   * Constructor that defines a string peptide with positions
   */
  SelfLoopPeptide(
    char* peptide, ///< character sequence of peptide
    int posA, ///< 1st link position on peptide
    int posB ///< 2nd link position on peptide
  );

  /**
   * Default destructor
   */
  virtual ~SelfLoopPeptide() {}

  /**
   * Adds Self-loop candidates to the collection
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass
    FLOAT_T max_mass, ///< max mass
    bool is_decoy,
    XLinkMatchCollection& candidates ///< collection to add candidates
    );

  /**
   *  \returns self-loop candidate
   */
  virtual XLINKMATCH_TYPE_T getCandidateType();
  
  /**
   * \returns the sequence string, marking the link sites with (X,Y)
   */
  virtual std::string getSequenceString();
  
  /**
   * \returns the mass of the self-loop peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< AVERAGE or MONO
  );

  /**
   * \returns a shuffled version of self-loop candidate
   */
  virtual void shuffle(
    std::vector<XLinkMatch*>& decoys
  );
  
  /**
   *  Predictes the ions for the self loop candidate
   */
  virtual void predictIons(
    IonSeries* ion_series, ///< ion vector to place ions in
    int charge ///< charge state
  );
  
  /**
   * \returns the sequence of the ion
   */
  std::string getIonSequence(
    Ion* ion ///< ion -in
  );

  /**
   *\returns the peptide object of self-loop peptide (0 is only valid)
   */
  virtual Crux::Peptide* getPeptide(
    int peptide_idx
  );

  /**
   *\returns the number of missed cleavages
   */
  virtual int getNumMissedCleavages();

  /**
   *\returns whether the peptide is modified by a variable mod
   */
  virtual bool isModified();
  
  
  virtual std::string getUnshuffledSequence();
  
  SelfLoopPeptide* getUnshuffledTarget();

  /**
   *\returns "target" or "decoy"
   */
  string getDecoyType();

};

bool compareSelfLoopPeptideMass(
				const SelfLoopPeptide& spep1,
				const SelfLoopPeptide& spep2);
bool compareSelfLoopPeptideMassToFLOAT(const SelfLoopPeptide& spep1, FLOAT_T mass);


#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
