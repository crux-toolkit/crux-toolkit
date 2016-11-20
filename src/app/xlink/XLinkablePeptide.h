/**
 * \file XLinkablePeptide.h
 * $Revision: 1.00 $
 * \brief Object for finding and defining the link sites on a peptide.
 *************************************************************************/
#ifndef XLINKABLEPEPTIDE_H_
#define XLINKABLEPEPTIDE_H_

#include "model/objects.h"
#include "model/Peptide.h"
#include "util/CacheableMass.h"
#include <vector>
#include <string>

#include "XLinkBondMap.h"

/**
 * \class XLinkablePeptide
 * \brief object for finding and defining the link sites on a peptide
 */
class XLinkablePeptide : public CacheableMass {
 protected:
  Crux::Peptide* peptide_; ///< the peptide object of this XLinkablePeptide (can be null)
  XLinkablePeptide* decoy_; ///< a saved decoy of the XLinkablePeptide.
  char* sequence_; ///< the sequence 
  std::vector<int> link_sites_; ///< the sequence indices where linking is possible
  bool is_decoy_; //Is this from the decoy database?
  size_t xcorr_link_idx_;
  FLOAT_T xcorr_;
  int index_;

  MODIFIED_AA_T* mod_seq_;

  int predict_ions_call_count_;
  /**
   * Initialize object
   */
  void init();
  
 public:
  
  /**
   * Default constructor
   */
  XLinkablePeptide();

  /**
   *  Constructor that defines the peptide as a sequence string
   */
  XLinkablePeptide(
    char* sequence ///< the peptide sequence
    );

  XLinkablePeptide(
    const XLinkablePeptide& xlinkablepeptide
  );

  XLinkablePeptide(
    XLinkablePeptide& xlinkablepeptide
  );

  /**
   * Constructor that defines the peptide and the linking sites
   */
  XLinkablePeptide(
    Crux::Peptide* peptide, ///< the peptide object 
    std::vector<int>& link_sites ///< the linking sites
    );

  /**
   * Constructor that generates the possible linking sites
   * on the peptide using the bondmap object
   */
  XLinkablePeptide(
    Crux::Peptide* peptide, ///< the peptide object 
    XLinkBondMap& bondmap, ///< the bond map
    int additional_cleavages = 0 ///< 1 - for selfloop peptides
  );
  
  /**
   * Default destructor
   */
  virtual ~XLinkablePeptide();
 
  /**
   * given a peptide and a XLinkBondMap object,
   * generate the possible link sites
   */
  static void findLinkSites(
    Crux::Peptide* peptide,  ///< the peptide object -in
    XLinkBondMap& bondmap,  ///< the bond map -in 
    std::vector<int>& link_sites, ///< the found link sites -out
    int additional_cleavages=0 ///< 0 for xlinks, 1 for selfloops
  );

  void clearSites();

  /**
   * \returns whether a link at this index in the sequence
   * would prevent cleavage
   */
  bool linkSeqPreventsCleavage(
    int seq_idx ///< the sequence index
  );

  static bool linkSeqPreventsCleavage(
    Crux::Peptide* peptide, ///< the peptide
    int seq_idx ///< the sequence index
  );

  /**
   * \returns whether the link site prevents enzymatic
   * cleavage
   */
  bool linkSitePreventsCleavage(
    size_t link_site_idx ///< link site
  );

  /**
   * \returns whether the peptide is cleavable at this index
   */
  static bool isCleavageSite(
    Crux::Peptide* peptide, ///< the peptide
    int seq_idx ///< the sequence index
  );

  /**
   * \returns whether the aa at this index is a missed cleavage
   */
  static bool isMissedCleavageSite(
    Crux::Peptide* peptide, ///< the peptide
    int seq_idx ///< the sequence index
  );

  /**
   * shuffle the linkable peptide, preserving the link sites
   * \returns the shuffled XlinkablePeptide
   */
  XLinkablePeptide shuffle();

  XLinkablePeptide* getCachedDecoy();
  

  /**
   * \returns the number of link sites on this peptide
   */
  size_t numLinkSites();

  /**
   * \returns whether this peptide is linkable or not
   */
  bool isLinkable();
  

  void setDecoy(bool is_decoy);

  /**
   * \returns whether the peptide is a decoy or not
   */
  bool isDecoy();


  /**
   * \returns whether the peptide is linkable using the bond map
   */
  static bool isLinkable(
    Crux::Peptide* peptide, ///< the peptide object
    XLinkBondMap& bondmap ///< the bond map
  );

  /**
   * Adds a link site by sequence index
   * \returns the index of the added link site
   */
  int addLinkSite(
    int seq_idx ///< the sequence index of the link site
  );

  /**
   * \retuns the sequence index of a link site
   */
  int getLinkSite(
    int link_site_idx ///< the index of the link site
    );

  /**
   * \returns the number of missed cleavages in the linkable peptide
   */
  int getMissedCleavageSites();
  static int getMissedCleavageSites(
    Crux::Peptide* peptide ///< the peptide
  );

  /**
   * \returns the peptide object associated with this XLinkablePeptide
   */
  Crux::Peptide* getPeptide();

  /**
   * \returns the mass of the xlinkable peptide
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type=MONO ///< MONO or AVERAGE
  );

  /**
   * \returns an allocated sequence c-string.  Must be freed
   */
  const char* getSequence();

  /**
   * \returns the allocated modified sequence. Must be freed.
   */
  MODIFIED_AA_T* getModifiedSequence();

  const MODIFIED_AA_T* getModifiedSequencePtr();
  /**
   * \returns the modified sequence string of the xlinkable peptide
   */
  std::string getModifiedSequenceString();

  
  void predictIons(
    IonSeries* ion_series,
    int charge,
    int link_idx,
    FLOAT_T mod_mass,
    bool clear = true
    );

  FLOAT_T getXCorr() const;  
  void setXCorr( 
    size_t link_idx,
    FLOAT_T xcorr
  );

  bool hasXCorr() const;
  
  
  /**
   * Is the linkable peptide modified?
   */
  bool isModified();
  /**
   * \returns whether the peptide is less than (by lexical modified sequence)
   */
  bool operator < (
    XLinkablePeptide other ///< the other XLinkablePeptide to compare to.
  ) const;

  void setIndex(int idx) {index_ = idx;}
  int getIndex() {return index_;}

  
};

bool compareXLinkablePeptideMass(const XLinkablePeptide& xpep1, const XLinkablePeptide& xpep2);
bool compareXLinkablePeptideMassToFLOAT(
  const XLinkablePeptide& xpep1, 
  const FLOAT_T& mass);

bool compareXLinkablePeptideMassToFLOAT2(
					 const FLOAT_T& mass,
					 const XLinkablePeptide& xpep1);


bool compareXLinkableXCorr(
  const XLinkablePeptide& xpep1,
  const XLinkablePeptide& xpep2
  );

bool compareXLinkableXCorrPtr(
  const XLinkablePeptide* xpep1,
  const XLinkablePeptide* xpep2
  );

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
