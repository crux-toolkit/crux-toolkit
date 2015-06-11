/**
 * \file SelfLoopPeptide.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a self-loop peptide in an xlink search
 *****************************************************************************/
#ifndef SELFLOOPPEPTIDE_H_
#define SELFLOOPPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"

#include <vector>

class SelfLoopPeptide : public XLinkMatch {
 protected:

  XLinkablePeptide linked_peptide_; ///< linkable peptide
  std::vector<int> link_pos_idx_; ///< the two link indices on the peptide
  bool is_decoy_; ///< indicate whether this is a decoy
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
  virtual ~SelfLoopPeptide() {};

  /**
   * Adds Self-loop candidates to the collection
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass
    FLOAT_T max_mass, ///< max mass
    XLinkBondMap& bondmap,  ///< valid link sites
    Index* index, ///< protein index
    Database* database, ///< protein database
    PEPTIDE_MOD_T** peptide_mods, ///< allowable modifications
    int num_peptide_mods, ///< number of allowable modifications
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
  virtual XLinkMatch* shuffle();
  
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

};

#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
