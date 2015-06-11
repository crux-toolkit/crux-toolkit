/**
 * \file XLinkMatchCollection.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Collection of possible xlink products
 *****************************************************************************/

#include "XLinkMatchCollection.h"
#include "XLinkPeptide.h"
#include "LinearPeptide.h"
#include "SelfLoopPeptide.h"
#include "XLinkScorer.h"

#include "Spectrum.h"

#include <iostream>


static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.0;       // For now, turn off the threshold.
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate, 
  int isotope,
  FLOAT_T window,
  WINDOW_TYPE_T precursor_window_type,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {

  //cerr <<"mz: "
  //     <<precursor_mz
  //     <<" charge:"
  //     <<charge
  //     <<" mass:"<<mass
  //     <<" window:"<<window<<endl;
  if (precursor_window_type == WINDOW_MASS) {
    //cerr<<"WINDOW_MASS"<<endl;
    min_mass = zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON - window;
    max_mass = zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    //cerr<<"WINDOW_MZ"<<endl;
    double min_mz = precursor_mz - window;
    double max_mz = precursor_mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge();
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge();
  } else if (precursor_window_type == WINDOW_PPM) {
    //cerr<<"WINDOW_PPM"<<endl;
    min_mass = (zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON) / (1.0 + window * 1e-6);
    max_mass = (zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON)/ (1.0 - window * 1e-6);
  }
  
  //cerr<<"min:"<<min_mass<<" "<<"max: "<<max_mass<<endl;

}

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  int isotope,
  bool use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass) {
  
  if (use_decoy_window) {
    get_min_max_mass(precursor_mz,
         zstate,
                     isotope,
         get_double_parameter("precursor-window-weibull"),
         get_window_type_parameter("precursor-window-type-weibull"),
         min_mass,
         max_mass);
  } else {
    get_min_max_mass(precursor_mz,
         zstate,
                     isotope,
         get_double_parameter("precursor-window"),
         get_window_type_parameter("precursor-window-type"),
         min_mass,
         max_mass);
  }
}

/**
 * Default constructor
 */
XLinkMatchCollection::XLinkMatchCollection() : MatchCollection() {
  carp(CARP_DEBUG, "XLinkMatchCollection():start");
  scan_ = 0;
}

/**
 * Copy constructor
 */
XLinkMatchCollection::XLinkMatchCollection(
  XLinkMatchCollection& vector ///<collection to copy
  ) : MatchCollection() {
  
  carp(CARP_DEBUG, "XLinkMatchCollection(XLinkMatchCollection):start");
  precursor_mz_ = vector.precursor_mz_;
  zstate_ = vector.zstate_;
  scan_ = vector.scan_;

  for (int idx=0;idx<vector.getMatchTotal();idx++) {
    XLinkMatch* currentCandidate = (XLinkMatch*)vector[idx];
    XLinkMatch* copyCandidate = NULL;
    switch (currentCandidate -> getCandidateType()) {
    case LINEAR_CANDIDATE:
    case DEADLINK_CANDIDATE:
      copyCandidate = 
  new LinearPeptide(*(LinearPeptide*)currentCandidate);
      break;
    case SELFLOOP_CANDIDATE:
      copyCandidate =
  new SelfLoopPeptide(*(SelfLoopPeptide*)currentCandidate);
      break;
    case XLINK_INTER_CANDIDATE:
    case XLINK_INTRA_CANDIDATE:
    case XLINK_INTER_INTRA_CANDIDATE:
      copyCandidate =
  new XLinkPeptide(*(XLinkPeptide*)currentCandidate);
      break;
    }
    add(copyCandidate);
  }


}

/**
 * Constructor that finds all possible candidates
 */
XLinkMatchCollection::XLinkMatchCollection(
  XLinkBondMap& bondmap, ///< allowable links
  PEPTIDE_MOD_T** peptide_mods, ///< list of possible modifications
  int num_peptide_mods, ///< number of possible modifications
  Index* index, ///< protein index
  Database* database ///< protein database
  ) {
  
  carp(CARP_DEBUG, "XLinkMatchCollection(...)");

  FLOAT_T min_mass = get_double_parameter("min-mass");
  FLOAT_T max_mass = get_double_parameter("max-mass");

  addCandidates(
    min_mass, 
    max_mass, 
    bondmap, 
    index, 
    database, 
    peptide_mods, 
    num_peptide_mods);
  
}

/**
 * Constructor that finds all candidates within a mass range
 */
void XLinkMatchCollection::addCandidates(
  FLOAT_T min_mass, ///< min mass
  FLOAT_T max_mass, ///< max mass
  XLinkBondMap& bondmap, ///< allowable links
  Index* index, ///< protein index
  Database* database, ///< protein database
  PEPTIDE_MOD_T** peptide_mods, ///< list of possible modifications
  int num_peptide_mods ///< number of possible modifications
  ) {

  carp(CARP_DEBUG, "XLinkMatchCollection.addCandidates() start");

  include_linear_peptides_ = get_boolean_parameter("xlink-include-linears");
  include_self_loops_ = get_boolean_parameter("xlink-include-selfloops");

  carp(CARP_DEBUG, "Adding xlink candidates");

  XLinkPeptide::addCandidates(
    min_mass, 
    max_mass,
    bondmap,
    index,
    database,
    peptide_mods,
    num_peptide_mods,
    *this);

  if (include_linear_peptides_) {

    LinearPeptide::addCandidates(
      min_mass,
      max_mass,
      index,
      database,
      peptide_mods,
      num_peptide_mods,
      *this);

  }

  if (include_self_loops_) {
  
    SelfLoopPeptide::addCandidates(
      min_mass,
      max_mass,
      bondmap,
      index,
      database,
      peptide_mods,
      num_peptide_mods,
      *this);
  }
}


/**
 * Constructor for finding all candidates within a mass range
 */
XLinkMatchCollection::XLinkMatchCollection(
  FLOAT_T precursor_mz,  ///< precursor m/z
  SpectrumZState& zstate, ///< z-state
  XLinkBondMap& bondmap, ///< allowable links
  Index* index,  ///< protein index
  Database* database,  ///< protein database
  PEPTIDE_MOD_T** peptide_mods,  ///< list of allowable peptide mods
  int num_peptide_mods,  ///< number of allowable peptides
  bool use_decoy_window  
  ) {

  carp(CARP_DEBUG, "Inside XLinkMatchCollection....");

  precursor_mz_ = precursor_mz;
  setZState(zstate);  


  FLOAT_T min_mass;
  FLOAT_T max_mass;
  vector<int> isotopes = get_int_vector_parameter("isotope-windows");
  for (int idx = 0; idx < isotopes.size();idx++) {
    get_min_max_mass(precursor_mz, zstate, isotopes[idx], use_decoy_window, min_mass, max_mass);
    carp(CARP_INFO, "isotope %i min:%g max:%g", isotopes[idx], min_mass, max_mass);
    addCandidates(min_mass, max_mass, bondmap, index, database, peptide_mods, num_peptide_mods);
  }
}

/**
 * adds a candidate to the list
 */
void XLinkMatchCollection::add(
  XLinkMatch* candidate ///< candidate to add
  ) {

  candidate->setZState(zstate_);
  candidate->setParent(this);
  addMatch(candidate);
  candidate->decrementPointerCount();
  experiment_size_++;

}

/**
 * \returns a candidate from the list by index
 */
XLinkMatch* XLinkMatchCollection::at(
  int idx ///< index of the candidate
  ) {
  if (idx < 0 || idx >= getMatchTotal()) {
    carp(CARP_FATAL, "XLinkMatchCollection:index %d out of bounds (0,%d)",
      idx, getMatchTotal());
  }
  return (XLinkMatch*)match_[idx];
}

/**
 * \returns a candidate from the list by index
 */
XLinkMatch* XLinkMatchCollection::operator [](
  int idx ///< index of the candidate
  ) {
  return (XLinkMatch*)match_[idx];
}


/**
 * shuffles the candidates and places the results in a decoy collection
 */
void XLinkMatchCollection::shuffle(
  XLinkMatchCollection& decoy_vector ///< collection to add decoys to
  ) {
  
  decoy_vector.precursor_mz_ = precursor_mz_;
  decoy_vector.zstate_ = zstate_;
  decoy_vector.scan_ = scan_;

  for (int idx=0;idx<getMatchTotal();idx++) {
    decoy_vector.add(at(idx)->shuffle());
  }

}

/**
 * scores all candidates against the spectrum
 */
void XLinkMatchCollection::scoreSpectrum(
  Crux::Spectrum* spectrum ///< spectrum to score against
  ) {

  int max_ion_charge = get_max_ion_charge_parameter("max-ion-charge");

  carp(CARP_DEBUG, "Creating scorer");
  XLinkScorer scorer(
    spectrum, 
    min(zstate_.getCharge(), max_ion_charge));

  for (int idx=0;idx<getMatchTotal();idx++) {
    carp(CARP_DEBUG, "Scoring candidate:%d", idx);
    scorer.scoreCandidate(at(idx));
  }

  // set the match_collection as having been scored
  scored_type_[XCORR] = true;
  if (get_boolean_parameter("compute-sp")) {
    scored_type_[SP] = true;
  }

  carp(CARP_DEBUG, "Done scoreSpectrum");
}

/**
 * fits a weibull to the collection
 */
void XLinkMatchCollection::fitWeibull() {

  //create the array of x's and 
  shift_=0;
  eta_=0;
  beta_=0;
  correlation_=0;

  FLOAT_T* xcorrs = extractScores(XCORR);
// reverse sort the scores

  std::sort(xcorrs, xcorrs + getMatchTotal(), greater<FLOAT_T>());

  double fraction_to_fit = get_double_parameter("fraction-top-scores-to-fit");
  int num_tail_samples = (int)(getMatchTotal() * fraction_to_fit);

  fit_three_parameter_weibull(xcorrs,
            num_tail_samples,
            getMatchTotal(),
            MIN_XCORR_SHIFT,
            MAX_XCORR_SHIFT,
            XCORR_SHIFT,
            CORR_THRESHOLD,
            &eta_,
            &beta_,
            &shift_,
            &correlation_);

  myfree(xcorrs);

}

/**
 * computes the p-value for the candidate
 */
void XLinkMatchCollection::computeWeibullPValue(
  int idx ///< candidate
  ) {

  at(idx)->computeWeibullPvalue(shift_, eta_, beta_);
}

/**
 * sets the scan for the collection
 */
void XLinkMatchCollection::setScan(
  unsigned int scan ///< scan number to set
  ) {
  scan_ = scan;
}

/**
 *\returns the scan number
 */
unsigned int XLinkMatchCollection::getScan() {
  return scan_;
}

/**
 *\returns the charge state of the collection
 */
int XLinkMatchCollection::getCharge() {
  return zstate_.getCharge();
}

/**
 *\returns the precursor m/z
 */
FLOAT_T XLinkMatchCollection::getPrecursorMZ() {
  return precursor_mz_;
}

/**
 * \returns the neutral mass of the collection
 */
FLOAT_T XLinkMatchCollection::getSpectrumNeutralMass() {
  return zstate_.getNeutralMass();
}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */

