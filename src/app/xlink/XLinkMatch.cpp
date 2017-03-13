/**
 * \file XLinkMatch.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a Match in an xlink search
 *****************************************************************************/
#include "XLinkMatch.h"

#include "parameter.h"
#include "model/Scorer.h"
#include "util/GlobalParams.h"
#include <sstream>
#include <ios>
#include <iomanip>
#include <iostream>
#include "XLinkPeptide.h"
#include "XLinkablePeptide.h"
#include "io/OutputFiles.h"
using namespace std;

/**
 * Constructor for XLinkMatch
 */
XLinkMatch::XLinkMatch() : Match() {
  cached_sequence_ = "";
  parent_ = NULL;
  pvalue_= 1;
  for (int idx = 0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = false;
    mass_[idx] = 0;
  }
  is_decoy_ = false;
}

/**
 * Default destrcutor for XLinkMatch
 */
XLinkMatch::~XLinkMatch() {
  //TODO deleted cached ions
}

void XLinkMatch::decrementPointerCount() {

  pointer_count_--;
}

/**
 * computes the pvalue for this match using the provided weibull paramters
 */
void XLinkMatch::computeWeibullPvalue(
  FLOAT_T shift, ///< shift parameter for weibull
  FLOAT_T eta, ///< eta parameter for weibull
  FLOAT_T beta ///< beta parameter for weibull
  ) {

  pvalue_ = compute_weibull_pvalue(getScore(XCORR), eta, beta, shift);
}

void XLinkMatch::setPValue(FLOAT_T pvalue) {
  pvalue_ = pvalue;
}

FLOAT_T XLinkMatch::getPValue() {

  return(pvalue_);

}

/**
 * \returns the protein id string for this match
 * default, can be overridden
 */
string XLinkMatch::getProteinIdString() {
  Crux::Peptide* peptide = this -> getPeptide(0);

  if (peptide == NULL) {
    return string("");
  } else {
    return XLink::get_protein_ids_locations(peptide);
  }
}

/**
 *\returns the protein id for this match
 */
string XLinkMatch::getProteinIdXString() {
  Crux::Peptide* peptide = this -> getPeptide(0);
  if (peptide == NULL) {
    return string("");
  } else {
    return XLink::get_protein_ids_locations(peptide);
  }

}

string XLinkMatch::getUnshuffledSequence() {
  Crux::Peptide* peptide = this->getPeptide(0);
  if (peptide == NULL) {
    return string("");
  } else {
    return peptide->getUnshuffledSequence();
  }
}


/**
 *\returns the flanking amino acids for the match
 */
string XLinkMatch::getFlankingAAString() {
  Crux::Peptide* peptide = this -> getPeptide(0);

  string ans("");

  if (peptide != NULL) {
    char* flanking_aas = peptide->getFlankingAAs();
    ans = flanking_aas;
    free(flanking_aas);
  }
  return ans;
}

/**
 * \returns the candidate type from a string value
 */
XLINKMATCH_TYPE_T XLinkMatch::getCandidateType(std::string& candidate) {
//  cerr << "get candidate type:"<<candidate<<endl;
  if (candidate == string("linear")) {
    return LINEAR_CANDIDATE;
  } else if (candidate == string("dead-link")) {
    return DEADLINK_CANDIDATE;
  } else if (candidate == string("self-loop")) {
    return SELFLOOP_CANDIDATE;
  } else if (candidate == string("xlink-inter")) {
    return XLINK_INTER_CANDIDATE;
  } else if (candidate == string("xlink-intra")) {
    return XLINK_INTRA_CANDIDATE;
  } else if (candidate == string("xlink-inter-intra")) {
    return XLINK_INTER_INTRA_CANDIDATE;
  } else {
    return INVALID_CANDIDATE;
  }

}

/**
 *\returns the string value for the candidate tyoe of the match
 */
string XLinkMatch::getCandidateTypeString() {
  return getCandidateTypeString(getCandidateType());
}

const string& XLinkMatch::getSequenceStringConst() {
  if (cached_sequence_ == "") {
    cached_sequence_ = getSequenceString();
  }
  return(cached_sequence_);
}


/**
 *\returns the string value of the given candidate type
 */
string XLinkMatch::getCandidateTypeString(
  XLINKMATCH_TYPE_T candidate ///< candidate type to convert
  ) {

  string ans = "";
  switch(candidate) {
    case LINEAR_CANDIDATE:
      ans = "linear";
      break;
    case DEADLINK_CANDIDATE:
      ans = "dead-link";
      break;
    case SELFLOOP_CANDIDATE:
      ans = "self-loop";
      break;
    case XLINK_INTER_CANDIDATE:
      ans = "xlink-inter";
      break;
    case XLINK_INTRA_CANDIDATE:
      ans = "xlink-intra";
      break;
    case XLINK_INTER_INTRA_CANDIDATE:
      ans = "xlink-inter-intra";
      break;
    default:
      ans = "unknown";
  }

  return ans;


}

vector<IonConstraint*> XLinkMatch::ion_constraint_xcorr_;

IonConstraint* XLinkMatch::getIonConstraintXCORR(int charge) {
  int idx = charge-1;
  while(ion_constraint_xcorr_.size() < charge) {
    ion_constraint_xcorr_.push_back(NULL);
  }
  if (ion_constraint_xcorr_[idx] == NULL) {
    ion_constraint_xcorr_[idx] = 
      IonConstraint::newIonConstraintSmart(XCORR, charge);
  } else {
    //    carp(CARP_INFO, "IonConstraint cache hit!");
  }
  return(ion_constraint_xcorr_[idx]);
}

IonSeries* XLinkMatch::getIonSeriesXCORR(int charge) {

  int idx = charge-1;
  while(ion_series_xcorr_.size() < charge) {
    ion_series_xcorr_.push_back(NULL);
  }

  if (ion_series_xcorr_[idx] == NULL) {
    IonSeries* ion_series_xcorr = new IonSeries(getIonConstraintXCORR(charge), charge); 
    predictIons(ion_series_xcorr, charge);
    ion_series_xcorr_[idx] = ion_series_xcorr;
  } else {
    carp(CARP_INFO, "Cache hit!");
  }
  
  return(ion_series_xcorr_[idx]);

}

/**
 * \returns the mass error in part-per-million (ppm)
 */
FLOAT_T XLinkMatch::getPPMError() {
  FLOAT_T mono_mass = getMass(MONO);
  FLOAT_T obs_mass = parent_->getSpectrumNeutralMass();
  FLOAT_T isotope;
  
  if (mono_mass > obs_mass) {
    isotope = floor((mono_mass - obs_mass) / MASS_NEUTRON + 0.5);
  } else {
    isotope = -floor((obs_mass - mono_mass) / MASS_NEUTRON + 0.5);
  }
  

  obs_mass = obs_mass + isotope * MASS_NEUTRON;

  FLOAT_T ppm = (mono_mass - obs_mass) / obs_mass * 1e6;
  return ppm;
  
}

/**
 * sets the XLinkMatchCollection owner of the match
 */
void XLinkMatch::setParent(XLinkMatchCollection* parent) {
  parent_ = parent;
}

const vector<XLinkMatch*>& XLinkMatch::getDecoys() {
  decoys_.clear();
  shuffle(decoys_);
  return(decoys_);
  /*
  if (decoys_.size() == 0) {
    shuffle(decoys_);
    if (decoys_.size() == 0) {
      carp(CARP_FATAL, "Can't get decoy(s) for XLinkMatch!");
    }
  } else {
    carp(CARP_INFO, "XLinkMatch::Using cached decoys");
  }
  
  return(decoys_);
  */
}


/**
 * Print one field in the tab-delimited output file, based on column index.
 */
void XLinkMatch::printOneMatchField(
  int      column_idx,             ///< Index of the column to print. -in
  MatchCollection* collection,  ///< collection holding this match -in 
  MatchFileWriter*    output_file,            ///< output stream -out
  Crux::Spectrum* spectrum, 
  int      num_target_matches,            ///< target matches for this spectrum -in
  int      num_decoy_matches, ///< decoy matches for this spectrum -in
  int      b_y_total,              ///< total b/y ions -in
  int      b_y_matched             ///< Number of b/y ions matched. -in
) {

  carp(CARP_DETAILED_DEBUG, "XLinkMatch::printOneMatchField:%s", get_column_header(column_idx));

  switch ((MATCH_COLUMNS_T)column_idx) {

  case PEPTIDE_MASS_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, getMass(GlobalParams::getIsotopicMass()));
    break;
  case PVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, pvalue_);
    break;
  case SEQUENCE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getSequenceStringConst());
    break;
  case PROTEIN_ID_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx, 
      getProteinIdString()); 
    break;
  case FLANKING_AA_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getFlankingAAString());
    break;
  case TARGET_DECOY_COL:
    switch (getCandidateType()) {
      case LINEAR_CANDIDATE:
      case DEADLINK_CANDIDATE:
      case SELFLOOP_CANDIDATE:
        if (this->getPeptide(0)->isDecoy()) {
          output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                           "decoy");
        } else {
          output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                           "target");
        }
        break;
      case XLINK_INTER_CANDIDATE:
      case XLINK_INTRA_CANDIDATE:
      case XLINK_INTER_INTRA_CANDIDATE:
        output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx,
                                         getDecoyType());
        break;
      case INVALID_CANDIDATE:
        carp(CARP_FATAL, "Encountered invalid candidate type.");
    }
    break;
  case XLINK_PRODUCT_TYPE_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getCandidateTypeString());
    break;
  case PPM_ERROR_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getPPMError());
    break;
  case XCORR_FIRST_COL:
    if ((Params::GetInt("xlink-top-n") != 0) &&
        (getCandidateType() == XLINK_INTER_CANDIDATE || 
        getCandidateType() == XLINK_INTRA_CANDIDATE || 
        getCandidateType() == XLINK_INTER_INTRA_CANDIDATE)) {
      XLinkPeptide *xpep = (XLinkPeptide*)this;
      XLinkablePeptide& lpep = xpep->getXLinkablePeptide(0);
      output_file->setColumnCurrentRow(
        (MATCH_COLUMNS_T)column_idx,
        lpep.getXCorr());
    } else {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 0);
    }
    break;
  case XCORR_SECOND_COL:
    if ((Params::GetInt("xlink-top-n") != 0) && 
        (getCandidateType() == XLINK_INTER_CANDIDATE ||
        getCandidateType() == XLINK_INTRA_CANDIDATE ||
        getCandidateType() == XLINK_INTER_INTRA_CANDIDATE)) {
      XLinkPeptide *xpep = (XLinkPeptide*)this;
      XLinkablePeptide& lpep = xpep->getXLinkablePeptide(1);
      output_file->setColumnCurrentRow(
        (MATCH_COLUMNS_T)column_idx,
        lpep.getXCorr());
    } else {
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 0);
    }
    break;
  case PROTEIN_ID_X_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getProteinIdXString());
    break;
  case XLINK_TYPE_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getCandidateTypeString());
    break;
  case ORIGINAL_TARGET_SEQUENCE_COL:
    if (null_peptide_ == true || Params::GetBool("concat")) {
      string seq = getUnshuffledSequence();
      output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, seq);
    }
    break;
  case MODIFICATIONS_COL:
    //TODO FIX!
    break;
  case ENZ_INT_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, getNumMissedCleavages());
    break;
  default:
    Match::printOneMatchField(column_idx,
      collection,
      output_file,
      spectrum, 
      num_target_matches,
      num_decoy_matches,
      b_y_total,
      b_y_matched
    );
  }
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

