/**
 * \file XLinkPeptide.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September 2014
 * \brief Object for Defining a crosslinked peptide in an xlink search
 *****************************************************************************/
#ifndef XLINKPEPTIDE_H_
#define XLINKPEPTIDE_H_

#include "model/objects.h"
#include "util/utils.h"

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
  
  XLinkPeptide* target_;
  
  static FLOAT_T pmin_;  ///< contains the minimum mass that a peptide from a crosslink product can take
  static bool pmin_set_; ///< has the pmin been set?
  /**
   * \returns the link position within each peptide
   */
  int getLinkPos(
    int peptide_idx ///< 0 - first peptide, 1 - second peptide
  );
 
  /*
   * Iterates through all linkable sites and adds valid xlink peptide candidates
   * \returns the number of candidates added.
   */
  static int addXLinkPeptides(
    XLinkablePeptide& pep1, ///< First linkable peptide
    XLinkablePeptide& pep2, ///< Second linkable peptide
    XLinkMatchCollection& candidates ///< XLinkable Candidates -out
  );

  /*
   * Is either of the participating peptides a decoy?
   */
  virtual bool isDecoy();

  /*
   * Distinguish between full targets, half target / half decoy, and
   * full decoy cross-linked pairs.
   */
  string getDecoyType();

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
  virtual ~XLinkPeptide() {}

  /**
   * makes sure that sequence1 is smaller in alphanumeric value than
   * sequence 2
   */
  void doSort();

  int getLinkIdx(
    int peptide_idx ///< 0 - first peptide, 1 -second peptide
  );

  
  
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
   * \returns whether the cross-link is from peptides that are both within the same protein and can be from different proteins
   */
  bool isInterIntra();



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
  static int addCandidates(
    FLOAT_T min_mass, ///< min mass of crosslinks
    FLOAT_T max_mass, ///< max mass of crosslinks
    vector<XLinkablePeptide>&, ///< 1st peptide iterator
    XLinkMatchCollection& candidates ///< candidates in/out
    );

  /**
   * adds crosslink candidates by iterating through all possible masses
   */
  static int addCandidates(
			    Crux::Spectrum* spectrum,
			    FLOAT_T precursor_mass,
			    int precursor_charge,
    FLOAT_T min_mass, ///< min mass of crosslinks
    FLOAT_T max_mass, ///< max mass of crosslinks
			    bool decoy,
    XLinkMatchCollection& candidates ///< candidates -in/out
    );

  /**
   * \returns the candidate type
   */
  virtual XLINKMATCH_TYPE_T getCandidateType();
  
  /**
   * \returns the sequence string
   */
  virtual std::string getSequenceString();
  
  virtual std::string getUnshuffledSequence();
  
  /**
   * \returns the mass of the xlink peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type ///< MONO or AVERAGE
  );

  /**
   * \returns a shuffled xlink peptide
   */
  virtual void shuffle(
    std::vector<XLinkMatch*>& decoys
  );

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

  XLinkablePeptide& getXLinkablePeptide(
    int peptide_idx ///<0 or 1
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
