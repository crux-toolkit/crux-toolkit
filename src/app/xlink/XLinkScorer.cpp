/**
 * \file XLinkScorer.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for scoring xlink products.  Tries to optimize the scoring time
 *****************************************************************************/

#include "XLinkScorer.h"
#include "model/Scorer.h"
#include "model/Ion.h"
#include "model/IonSeries.h"
#include "util/Params.h"
#include "XLinkPeptide.h"
#include "util/GlobalParams.h"


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

IonConstraint* XLinkScorer::getIonConstraintXCorr() {
  return ion_constraint_xcorr_;
}




/**
 * default constructor
 */
XLinkScorer::XLinkScorer() {

  init(NULL, 0, Params::GetBool("compute-sp"));
}

/**
 * Constructor for spectrum and charge
 */
XLinkScorer::XLinkScorer(
  Crux::Spectrum* spectrum, ///< spectrum to score
  int charge ///< charge state
  ) {

  init(spectrum, charge, Params::GetBool("compute-sp"));
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

  if (compute_sp_) {

    candidate->predictIons(ion_series_sp_, charge_);
    FLOAT_T sp = scorer_sp_->scoreSpectrumVIonSeries(spectrum_, ion_series_sp_);
    
    candidate->setScore(SP, sp);
    candidate->setScore(BY_IONS_MATCHED, scorer_sp_->getSpBYIonMatched());
    candidate->setScore(BY_IONS_TOTAL, scorer_sp_->getSpBYIonPossible());

  }
  FLOAT_T xcorr = 0;
  if (GlobalParams::getXLinkTopN() != 0 &&
      (candidate->getCandidateType() == XLINK_INTER_CANDIDATE ||
       candidate->getCandidateType() == XLINK_INTRA_CANDIDATE ||
       candidate->getCandidateType() == XLINK_INTER_INTRA_CANDIDATE
       )) {
        
    XLinkPeptide* xpep = (XLinkPeptide*)candidate;
    XLinkablePeptide& xpep1 = xpep->getXLinkablePeptide(0);
    XLinkablePeptide& xpep2 = xpep->getXLinkablePeptide(1);
    FLOAT_T link_mass = XLinkPeptide::getLinkerMass();
    FLOAT_T mod_mass1 = xpep1.getMass(GlobalParams::getIsotopicMass()) + link_mass;
    FLOAT_T mod_mass2 = xpep2.getMass(GlobalParams::getIsotopicMass()) + link_mass;
    size_t link_idx1 = xpep->getLinkIdx(0);
    size_t link_idx2 = xpep->getLinkIdx(1);
    
    FLOAT_T xcorr1;
    if (!xpep1.hasXCorr()) {    
      xcorr1 = scoreXLinkablePeptide(xpep1, link_idx1, mod_mass2);
      xpep1.setXCorr(link_idx1, xcorr1);
    } else {
      xcorr1 = xpep1.getXCorr();
    }
    
    FLOAT_T xcorr2;
    if (!xpep2.hasXCorr()) {
      xcorr2 = scoreXLinkablePeptide(xpep2, link_idx2, mod_mass1);
      xpep2.setXCorr(link_idx2, xcorr2);
    } else {
      xcorr2 = xpep2.getXCorr();
    }
        
    xcorr = xcorr1+xcorr2;
    candidate->setScore(XCORR, xcorr);
  } else {
    candidate->predictIons(ion_series_xcorr_, charge_);
    xcorr = scorer_xcorr_->scoreSpectrumVIonSeries(spectrum_, ion_series_xcorr_);
    candidate->setScore(XCORR, xcorr);
  }
  
  
/*
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
  */
  return xcorr;
}

FLOAT_T XLinkScorer::scoreXLinkablePeptide(
  XLinkablePeptide& xlpeptide,
  int link_idx,
  FLOAT_T mod_mass) {

  xlpeptide.predictIons(ion_series_xcorr_, charge_, link_idx, mod_mass);
  FLOAT_T xcorr = scorer_xcorr_->scoreSpectrumVIonSeries(spectrum_, ion_series_xcorr_);
  return xcorr;

}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
