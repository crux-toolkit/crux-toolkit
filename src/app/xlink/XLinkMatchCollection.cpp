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
#include "MonoLinkPeptide.h"
#include "XLinkScorer.h"

#include "model/Spectrum.h"
#include "util/GlobalParams.h"
#include "util/Params.h"
#include "util/StringUtils.h"

#include <iostream>


static const FLOAT_T MIN_XCORR_SHIFT = -5.0;
static const FLOAT_T MAX_XCORR_SHIFT  = 5.0;
//#define CORR_THRESHOLD 0.995   // Must achieve this correlation, else punt.
static const FLOAT_T CORR_THRESHOLD = 0.5;
static const FLOAT_T XCORR_SHIFT = 0.05;

using namespace std;

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate, 
  int isotope,
  FLOAT_T window,
  WINDOW_TYPE_T precursor_window_type,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass,
  FLOAT_T& precursor_mass) {

  //cerr <<"mz: "
  //     <<precursor_mz
  //     <<" charge:"
  //     <<charge
  //     <<" mass:"<<mass
  //     <<" window:"<<window<<endl;

  precursor_mass = zstate.getNeutralMass() + (double)isotope*MASS_NEUTRON;

  if (precursor_window_type == WINDOW_MASS) {
    //cerr<<"WINDOW_MASS"<<endl;
    min_mass = precursor_mass - window;
    max_mass = precursor_mass + window;
  } else if (precursor_window_type == WINDOW_MZ) {
    //cerr<<"WINDOW_MZ"<<endl;
    double min_mz = precursor_mz - window;
    double max_mz = precursor_mz + window;
    min_mass = (min_mz - MASS_PROTON) * (double)zstate.getCharge();
    max_mass = (max_mz - MASS_PROTON) * (double)zstate.getCharge();
  } else if (precursor_window_type == WINDOW_PPM) {
    //cerr<<"WINDOW_PPM"<<endl;
    min_mass = precursor_mass * (1.0 - window * 1e-6);
    max_mass = precursor_mass * (1.0 + window * 1e-6);
  }
  
  //cerr<<"min:"<<min_mass<<" "<<"max: "<<max_mass<<endl;

}

void get_min_max_mass(
  FLOAT_T precursor_mz, 
  SpectrumZState& zstate,
  int isotope,
  bool use_decoy_window,
  FLOAT_T& min_mass, 
  FLOAT_T& max_mass,
  FLOAT_T& precursor_mass) {
  
  if (use_decoy_window) {
    get_min_max_mass(precursor_mz,
         zstate,
                     isotope,
         Params::GetDouble("precursor-window-weibull"),
         string_to_window_type(Params::GetString("precursor-window-type-weibull")),
         min_mass,
		     max_mass, precursor_mass);
  } else {
    get_min_max_mass(precursor_mz,
         zstate,
                     isotope,
      GlobalParams::getPrecursorWindow(),
      GlobalParams::getPrecursorWindowType(),		     
      min_mass,
      max_mass, 
      precursor_mass);
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

  for (int idx = 0; idx < vector.getMatchTotal(); idx++) {
    XLinkMatch* currentCandidate = (XLinkMatch*)vector[idx];
    XLinkMatch* copyCandidate = NULL;
    switch (currentCandidate -> getCandidateType()) {
    case LINEAR_CANDIDATE:
      copyCandidate = 
      new LinearPeptide(*(LinearPeptide*)currentCandidate);
      break;
    case DEADLINK_CANDIDATE:
      copyCandidate = 
      new MonoLinkPeptide(*(MonoLinkPeptide*)currentCandidate);
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
    case INVALID_CANDIDATE:
      carp(CARP_ERROR, "Invalid candidate type.");
      exit(1);
    }
    add(copyCandidate);
  }


}

/**
 * Constructor that finds all possible candidates
 */
/*
XLinkMatchCollection::XLinkMatchCollection() {
  
  carp(CARP_DEBUG, "XLinkMatchCollection(...)");

  FLOAT_T min_mass = Params::GetDouble("min-mass");
  FLOAT_T max_mass = Params::GetDouble("max-mass");

  addCandidates(
		NULL,
		0,
		1,
		min_mass, 
		max_mass,
		false);

}
*/

/**
 * Constructor that finds all candidates within a mass range
 */
void XLinkMatchCollection::addCandidates(
  Crux::Spectrum *spectrum,
  FLOAT_T precursor_mass,
  int precursor_charge,
  FLOAT_T min_mass, ///< min mass
  FLOAT_T max_mass, ///< max mass
  bool decoy
  ) {

  carp(CARP_DETAILED_DEBUG, "XLinkMatchCollection.addCandidates() start");

  include_linear_peptides_ = GlobalParams::getXLinkIncludeLinears();
  include_self_loops_ = GlobalParams::getXLinkIncludeSelfloops();
  bool include_monolink_peptides = GlobalParams::getXLinkIncludeDeadends();
  if (GlobalParams::getXLinkIncludeInter() ||
      GlobalParams::getXLinkIncludeIntra() ||
      GlobalParams::getXLinkIncludeInterIntra()) {
    int num_xlink_candidates = 0;
    carp(CARP_DEBUG, "Adding xlink candidates");
    carp(CARP_DEBUG, "precursor:%g", precursor_mass);
    carp(CARP_DEBUG, "min:%g", min_mass);
    carp(CARP_DEBUG, "max:%g", max_mass);
    num_xlink_candidates = XLinkPeptide::addCandidates(
      spectrum,
      precursor_mass,
      precursor_charge,
      min_mass, 
      max_mass,
      decoy,
      *this);
    carp(CARP_DETAILED_DEBUG,"Number of xlink candidates:%d", num_xlink_candidates);
    if (num_xlink_candidates == 0 && Params::GetBool("require-xlink-candidate")) {
      carp(CARP_DEBUG, "no xlink candidate, returning");
      return;
    }
    
  }
  if (include_linear_peptides_) {
    carp(CARP_DETAILED_DEBUG, "adding Linear Candidates");
    LinearPeptide::addCandidates(
      min_mass,
      max_mass,
      decoy,
      *this);

  }
  
  if (include_monolink_peptides) {
    carp(CARP_DETAILED_DEBUG, "adding monolink candidates");
    MonoLinkPeptide::addCandidates(
      min_mass, max_mass, decoy, *this);
  }
  

  if (include_self_loops_) {
    carp(CARP_DETAILED_DEBUG, "adding selfloop candidates");
    SelfLoopPeptide::addCandidates(
      min_mass,
      max_mass,
      decoy,
      *this);
  }
  carp(CARP_DETAILED_DEBUG, "XLinkMatchCollection.addCandidates() done.");

}


/**
 * Constructor for finding all candidates within a mass range
 */
XLinkMatchCollection::XLinkMatchCollection(
  Crux::Spectrum *spectrum, ///< spectrum
  SpectrumZState& zstate, ///< z-state
  bool decoy,
  bool use_decoy_window  
  ) {

  carp(CARP_DEBUG, "Inside XLinkMatchCollection....");

  precursor_mz_ = spectrum->getPrecursorMz();
  setZState(zstate);  


  FLOAT_T min_mass;
  FLOAT_T max_mass;
  const vector<int>& isotopes = GlobalParams::getIsotopeWindows(); 
  for (int idx = 0; idx < isotopes.size();idx++) {
    FLOAT_T precursor_mass;
    get_min_max_mass(precursor_mz_, zstate, isotopes[idx], use_decoy_window, min_mass, max_mass, precursor_mass);
    carp(CARP_DEBUG, "isotope %i precursor: %g min:%g max:%g", isotopes[idx], precursor_mass, min_mass, max_mass);
    addCandidates(spectrum, precursor_mass, zstate.getCharge(), 
		  min_mass, max_mass, decoy);
  }
}

/**
 * adds a candidate to the list
 */
void XLinkMatchCollection::add(
  XLinkMatch* candidate, ///< candidate to add
  bool copy
  ) {

  candidate->setZState(zstate_);
  candidate->setParent(this);
  addMatch(candidate);
  if (!copy) {
    candidate->decrementPointerCount();
  }
  experiment_size_++;

}

void XLinkMatchCollection::add(
  const vector<XLinkMatch*>& candidates,
  bool copy
  ) {

  for (size_t idx = 0;idx < candidates.size();idx++) {
    //carp(CARP_INFO, "Adding %s", candidates.at(idx)->getSequenceString().c_str());
    add(candidates.at(idx), copy);
  }
  
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
XLinkMatch* XLinkMatchCollection::operator[] (
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

  for (int idx = 0; idx < getMatchTotal(); idx++) {
    decoy_vector.add(at(idx)->getDecoys());
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

  for (int idx = 0; idx < getMatchTotal(); idx++) {
    carp(CARP_DEBUG, "Scoring candidate:%d", idx);
    scorer.scoreCandidate(at(idx));
  }

  // set the match_collection as having been scored
  scored_type_[XCORR] = true;
  if (Params::GetBool("compute-sp")) {
    scored_type_[SP] = true;
  }

  carp(CARP_DEBUG, "Done scoreSpectrum");
}

/**
 * fits a weibull to the collection
 */
void XLinkMatchCollection::fitWeibull() {

  //create the array of x's and 
  shift_ = 0;
  eta_ = 0;
  beta_ = 0;
  correlation_ = 0;

  vector<FLOAT_T> xcorrs = extractScores(XCORR);

  // reverse sort the scores
  std::sort(xcorrs.begin(), xcorrs.end(), greater<FLOAT_T>());

  double fraction_to_fit = Params::GetDouble("fraction-top-scores-to-fit");
  int num_tail_samples = (int)(getMatchTotal() * fraction_to_fit);

  FLOAT_T* xcorr_array = new FLOAT_T[xcorrs.size()];
  std::copy(xcorrs.begin(), xcorrs.end(), xcorr_array);
  fit_three_parameter_weibull(xcorr_array,
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
  delete[] xcorr_array;
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

