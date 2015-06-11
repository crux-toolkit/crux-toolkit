/**
 * \file XLinkPeptide.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September 2014
 * \brief Object for Defining a crosslinked peptide in an xlink search
 *****************************************************************************/
#ifndef XLINKPEPTIDE_H_
#define XLINKPEPTIDE_H_

#include "objects.h"
#include "utils.h"

#include "XLinkMatch.h"
#include "XLinkBondMap.h"
#include "XLinkablePeptide.h"


#include <set>
#include <vector>

class XLinkPeptide : public XLinkMatch {
 protected:
  static FLOAT_T linker_mass_; ///< contains the link mass
  static std::set<Crux::Peptide*> allocated_peptides_; ///< tracks the peptides that have been allocated during the candidate finding step
  std::vector<XLinkablePeptide> linked_peptides_; ///< contains the two peptides in the match
  std::vector<int> link_pos_idx_; ///< contains the link positions for the two peptides
  
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///< indicator determining whether the mass has been calculated already
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///< hold the mass of the xlink peptide
  
  bool is_decoy_; ///< indicates whether the peptide is a decoy
  static FLOAT_T pmin_;  ///< contains the minimum mass that a peptide from a crosslink product can take
  static bool pmin_set_; ///< has the pmin been set?
  /**
   * \returns the link position within each peptide
   */
  int getLinkPos(
    int peptide_idx ///< 0 - first peptide, 1 - second peptide
  );

 public:
  
  /**
   * Default Constructor
   */
  XLinkPeptide();
  
  /**
   * Constructor using two linkable peptide objects and locations
   */
  XLinkPeptide(
    XLinkablePeptide& peptideA, ///< 1st peptide
    XLinkablePeptide& peptideB, ///< 2nd peptide
    int posA, ///< link pos1
    int posB  ///<link pos2
    );
  
  /**
   * Constructor for using two string peptides and locations
   */
  XLinkPeptide(
    char* peptideA, ///< sequence of peptide1
    char* peptideB, ///< sequence of peptide2
    int posA, ///< position of crosslink in peptide1
    int posB ///< position of crosslink in peptide2
    );

  /**
   * Destructor
   */
  virtual ~XLinkPeptide() {};

  /**
   * makes sure that sequence1 is smaller in alphanumeric value than
   * sequence 2
   */
  void doSort();

  /**
   * \returns whether the cross-link is from peptides from two different
   * proteins
   */
  bool isInter();
  
  /**
   * \returns whether the cross-link is from peptides within the same protein
   */
  bool isIntra();

  /**
   * sets the static linker mass variable
   */
  static void setLinkerMass(
    FLOAT_T linker_mass ///< linker mass
  );
  
  /**
   * \returns the linker mass from the static variable
   */
  static FLOAT_T getLinkerMass();

  /**
   * adds crosslink candidates to the XLinkMatchCollection using
   * the passed in iterator for the 1st peptide
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass of crosslinks
    FLOAT_T max_mass, ///< max mass of crosslinks
    XLinkBondMap& bondmap, ///< valid crosslink map
    Index* index, ///< protein index
    Database* database, ///< protein database
    PEPTIDE_MOD_T* peptide_mod2, ///< modifications for the second peptide
    bool decoy2, ///< is the second peptide a decoy?
    XLinkablePeptideIterator& iter1, ///< 1st peptide iterator
    XLinkMatchCollection& candidates ///< candidates in/out
    );

  /**
   * adds crosslink candidates by iterating through all possible masses
   */
  static void addCandidates(
    FLOAT_T min_mass, ///< min mass of crosslinks
    FLOAT_T max_mass, ///< max mass of crosslinks
    XLinkBondMap& bondmap, ///< valid crosslink map
    Index* index,  ///< protein index
    Database* database, ///< protein database
    PEPTIDE_MOD_T** peptide_mods, ///< available variable mods
    int num_peptide_mods, ///< number of available modifications
    XLinkMatchCollection& candidates ///< candidates -in/out
    );

  /**
   * Gets all peptides that are linkable, i.e. have link sites
   */
  static void addLinkablePeptides(
    double min_mass, ///< min mass of peptides
    double max_mass, ///< max mass of peptides
    Index* index, ///< protein index
    Database* database, ///< protein database 
    PEPTIDE_MOD_T* peptide_mod, ///< modifications
    bool is_decoy, ///< are the peptides decoys
    XLinkBondMap& bondmap, ///<valid crosslink map
    std::vector<XLinkablePeptide>& linkable_peptides ///< list of linkable peptides -out
    );

  /**
   * \returns the candidate type
   */
  virtual XLINKMATCH_TYPE_T getCandidateType();
  
  /**
   * \returns the sequence string
   */
  virtual std::string getSequenceString();
  
  /**
   * \returns the mass of the xlink peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
  );

  /**
   * \returns a shuffled xlink peptide
   */
  virtual XLinkMatch* shuffle();

  /**
   * fills the ion series with the predicted ions for the cross linked candidate
   */ 
  virtual void predictIons(
    IonSeries* ion_series, ///< IonSeries object to fill
    int charge ///< charge state of the candidate 
  );

  /**
   * predicts the ions for xlinked peptide
   */
  void predictIons(
    IonSeries* ion_series, ///< IonSeries to fill
    int charge, ///< charge state of the peptide
    bool first ///< is this the first peptide?
    );
  
  /**
  * \returns the sequence from the ion
  */
  std::string getIonSequence(
    Ion* ion ///< pointer to the ion object
  );
  
  /**
   * \returns the Peptide object for the xlinked peptide
   */
  virtual Crux::Peptide* getPeptide(
    int peptide_idx ///< 0 or 1
  );

  /**
   * \returns the number of missed cleavages for the cross-linked peptide
   */
  virtual int getNumMissedCleavages();

  /**
   *\returns whether the cross-linked peptide is modified
   */
  virtual bool isModified();

  /**
   *\returns the protein id string for the xlinked peptide
   */
  virtual std::string getProteinIdString();
  
  /**
   * \returns the protein id strings where the (X) is the position
   * within the protein
   */
  virtual std::string getProteinIdsXLocations(
  int idx ///< peptide index (0 or 1)
  );
  
  /**
   * \returns the protein id string where the (X)s are the positions in
   * the proteins
   */
  virtual std::string getProteinIdXString();
  
  /**
   * \returns the flanking amino acids for both peptides, separated by ;
   */
  virtual std::string getFlankingAAString();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
