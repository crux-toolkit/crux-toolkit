/**
 * \file XLinkMatchCollection.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Collection of possible xlink products
 *****************************************************************************/
#ifndef XLINKMATCHCOLLECTION_H_
#define XLINKMATCHCOLLECTION_H_

/* Crux Includes */
#include "model/objects.h"
#include "model/MatchCollection.h"
#include "model/Database.h"
#include "model/Spectrum.h"
#include "util/modifications.h"
#include "model/SpectrumZState.h"

/* XLink Includes */
#include "XLink.h"
#include "XLinkMatch.h"

class XLinkMatchCollection : public MatchCollection {
 protected:

  bool include_linear_peptides_; ///< Include linear peptides 
  bool include_self_loops_; ///< include self loops
  int scan_; ///< scan number of the collection
  FLOAT_T precursor_mz_; ///< precursor m/z
  Crux::Spectrum* spectrum_; ///< spectrum object

  /**
   * Adds all of the possible candidates given the mass range
   */
  void addCandidates(
    Crux::Spectrum *spectrum, ///<spectrum
    FLOAT_T precursor_mass, ///< precursor mass
    int precursor_charge,
    FLOAT_T min_mass, ///< minimum mass
    FLOAT_T max_mass,  ///< maximum mass
    bool decoy ///< decoys?
    );

 public:

  /**
   * Default constructor
   */
  XLinkMatchCollection();

  /**
   * Copy constructor
   */
  XLinkMatchCollection(
    XLinkMatchCollection& vector ///< match collection to copy
  );

  /**
   * Constructor for finding all candidates within a mass range
   */
  XLinkMatchCollection(
    Crux::Spectrum* spectrum, ///< Spectrum
    SpectrumZState& zstate, ///< z-state
    bool decoy,
    bool use_decoy_window = false ///< decoys?
  );

  /**
   * Default destructor
   */
  virtual ~XLinkMatchCollection() {}

  /**
   * adds a candidate to the list
   */
  void add(XLinkMatch* candidate, bool copy = false);
  
  void add(
    const vector<XLinkMatch*>& candidates,
    bool copy = false
  );
  
  /**
   *\returns a candidate from the list by index
   */
  XLinkMatch* operator[] (
    int idx ///< index
  );
  
  /**
   *\returns a candidate from the list by index
   */
  XLinkMatch* at(
    int idx ///< index
  );
  
  /**
   * shuffles the candidates and places the results in a decoy collection
   */
  void shuffle(
    XLinkMatchCollection& decoy_vector ///< collection to add decoys to
  );

  /**
   * scores all candidates against the spectrum
   */
  void scoreSpectrum(
    Crux::Spectrum* spectrum ///< spectrum to score against
  );
  
  /**
   * sets the ranks for the candidates
   */
  void setRanks();
  
  /**
   * fits a weibull to the xcorrs in the collection
   */
  void fitWeibull();
  
  FLOAT_T getEta() {return eta_;}
  FLOAT_T getBeta() {return beta_;}
  FLOAT_T getShift() {return shift_;}
  FLOAT_T getCorrelation() {return correlation_;}
  

  /**
   * computes the p-values for the candidate
   */
  void computeWeibullPValue(
    int idx ///< candidate index
    );

  /**
   * sets the scan for the collection
   */
  void setScan(
    unsigned int scan ///< scan number to set
  );
  
  /**
   * \returns the scan number
   */
  unsigned int getScan();
  
  /**
   * \returns the charge state of the collection
   */
  int getCharge();
  
  /**
   * \returns the precursor m/z
   */
  FLOAT_T getPrecursorMZ();
  
  /**
   * \returns the neutral mass of the collection
   */
  FLOAT_T getSpectrumNeutralMass();

};

#endif

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
