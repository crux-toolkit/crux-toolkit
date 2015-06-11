/**
 * \file XLinkMatchCollection.h
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Collection of possible xlink products
 *****************************************************************************/
#ifndef XLINKMATCHCOLLECTION_H_
#define XLINKMATCHCOLLECTION_H_

/* Crux Includes */
#include "objects.h"
#include "MatchCollection.h"
#include "Index.h"
#include "Database.h"
#include "modifications.h"
#include "SpectrumZState.h"

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
    FLOAT_T min_mass, ///< minimum mass
    FLOAT_T max_mass, ///< maximum mass
    XLinkBondMap& bondmap, ///< map of valid links
    Index* index, ///< protein index
    Database* database, ///< protein database
    PEPTIDE_MOD_T** peptide_mods, ///< list of possible mods
    int num_peptide_mods ///< number of possible mods
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
   * Constructor that finds all possible candidates
   */
  XLinkMatchCollection(
    XLinkBondMap& bondmap, ///< allowable links
    PEPTIDE_MOD_T** peptide_mods, ///< list of allowable peptide mods
    int num_peptide_mods, ///< number of allowable peptide mods
    Index* index, ///< protein index
    Database* database ///< protein database
    );

  /**
   * Constructor for finding all candidates within a mass range
   */
  XLinkMatchCollection(
    FLOAT_T precursor_mz, ///< precursor m/z
    SpectrumZState& zstate, ///< z-state
    XLinkBondMap& bondmap, ///< allowable links
    Index* index, ///< protein index
    Database* database, ///protein database
    PEPTIDE_MOD_T** peptide_mods, ///< list of allowable peptide mods
    int num_peptide_mods, ///< number of allowable peptides
    bool use_decoy_window=false ///< decoys?
  );

  /**
   * Default destructor
   */
  virtual ~XLinkMatchCollection() {};

  /**
   * adds a candidate to the list
   */
  void add(XLinkMatch* candidate);

  /**
   *\returns a candidate from the list by index
   */
  XLinkMatch* operator [](
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
