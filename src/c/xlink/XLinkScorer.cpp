/**
 * \file XLinkScorer.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for scoring xlink products.  Tries to optimize the scoring time
 *****************************************************************************/

#include "XLinkScorer.h"
#include "Scorer.h"
#include "Ion.h"
#include "IonSeries.h"
#include "XLinkPeptide.h"

#include <iostream>

using namespace std;

/**
 * initializes the object with the spectrum
 * charge and appropriate objects
 */
void XLinkScorer::init(
  Crux::Spectrum* spectrum, ///< spectrum to score
  int charge, ///< charge state
  bool compute_sp ///< are we scoring sp?
  ) {

  scorer_xcorr_ = new Scorer(XCORR);
  scorer_sp_ = new Scorer(SP);
  spectrum_ = spectrum;
  charge_ = charge;
  compute_sp_ = compute_sp;

  if ((spectrum_ != NULL) && (charge_ > 0)) {

    ion_constraint_xcorr_ = 
      IonConstraint::newIonConstraintSmart(XCORR, charge_);

    ion_constraint_sp_ =
      IonConstraint::newIonConstraintSmart(SP, charge_);

    ion_series_xcorr_ = 
      new IonSeries(ion_constraint_xcorr_, charge_);

    ion_series_sp_ =
      new IonSeries(ion_constraint_sp_, charge_);
  } else {

    ion_constraint_xcorr_ = NULL;
    ion_constraint_sp_ = NULL;
    ion_series_xcorr_ = NULL;
    ion_series_sp_ = NULL;
  }
}

/**
 * default constructor
 */
XLinkScorer::XLinkScorer() {

  init(NULL, 0, get_boolean_parameter("compute-sp"));
}

/**
 * Constructor for spectrum and charge
 */
XLinkScorer::XLinkScorer(
  Crux::Spectrum* spectrum, ///< spectrum to score
  int charge ///< charge state
  ) {

  init(spectrum, charge, get_boolean_parameter("compute-sp"));
}

/**
 * Constructor for spectrum, charge, and sp
 */
XLinkScorer::XLinkScorer(
  Crux::Spectrum* spectrum, ///<spectrum to score
  int charge,  ///< charge state
  bool compute_sp ///< computing sp as well?
  ) {

  init(spectrum, charge, compute_sp);
}

/**
 * Destructor
 */
XLinkScorer::~XLinkScorer() {
  //carp(CARP_INFO, "XLinkScorer::~XLinkScorer()");
  delete ion_series_xcorr_;
  delete ion_series_sp_;
  delete ion_constraint_xcorr_;
  delete ion_constraint_sp_;
  delete scorer_xcorr_;
  delete scorer_sp_;
  
}

/**
 * \returns the xcorr score for the candidate and sets the sp if requested
 */
FLOAT_T XLinkScorer::scoreCandidate(
  XLinkMatch* candidate ///< candidate to score
  ) {

  candidate->predictIons(ion_series_xcorr_, charge_);
  FLOAT_T xcorr = scorer_xcorr_->scoreSpectrumVIonSeries(spectrum_, ion_series_xcorr_);
  candidate->setScore(XCORR, xcorr);

  if (compute_sp_) {

    candidate->predictIons(ion_series_sp_, charge_);
    FLOAT_T sp = scorer_sp_->scoreSpectrumVIonSeries(spectrum_, ion_series_sp_);
    
    candidate->setScore(SP, sp);
    candidate->setBYIonInfo(scorer_sp_);

  }

  if (candidate->getCandidateType() == XLINK_INTER_CANDIDATE || candidate->getCandidateType() == XLINK_INTRA_CANDIDATE) {
    carp(CARP_DEBUG, "Scoring xlink peptides individually");
    XLinkPeptide* xlink_match = (XLinkPeptide*)candidate;
    carp(CARP_DEBUG, "ions first");
    xlink_match->predictIons(ion_series_xcorr_, charge_, true);
    carp(CARP_DEBUG, "scoring first");
    FLOAT_T xcorr1 = scorer_xcorr_->scoreSpectrumVIonSeries(spectrum_, ion_series_xcorr_);
    carp(CARP_DEBUG, "first:%f", xcorr1);
    xlink_match->predictIons(ion_series_xcorr_, charge_, false);
    carp(CARP_DEBUG, "scoring second");
    FLOAT_T xcorr2 = scorer_xcorr_->scoreSpectrumVIonSeries(spectrum_, ion_series_xcorr_); 
    carp(CARP_DEBUG, "second:%f", xcorr2); 
    candidate->setScore(XCORR_FIRST, xcorr1);
    candidate->setScore(XCORR_SECOND, xcorr2);
  } else {

    candidate->setScore(XCORR_FIRST, xcorr);
    candidate->setScore(XCORR_SECOND, xcorr);
  }

  return xcorr;

}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
