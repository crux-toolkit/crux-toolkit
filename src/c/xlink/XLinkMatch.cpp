/**
 * \file XLinkMatch.h 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a Match in an xlink search
 *****************************************************************************/
#include "XLinkMatch.h"

#include "parameter.h"
#include "Scorer.h"
#include <sstream>
#include <ios>
#include <iomanip>
#include <iostream>


using namespace std;

/**
 * Constructor for XLinkMatch
 */
XLinkMatch::XLinkMatch() {
  parent_ = NULL;
  pvalue_= 1;
  for (int idx = 0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = false;
    mass_[idx] = 0;
  }
}

/**
 * Default destrcutor for XLinkMatch
 */
XLinkMatch::~XLinkMatch() {

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
  cerr <<"eta:"<<eta<<" beta:"<<beta<<" shift:"<<shift<<endl;
  cerr <<"xcorr:"<<getScore(XCORR)<<" pvalue:"<<pvalue_<<endl;
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

/**
 * \returns the mass of the match
 */
FLOAT_T XLinkMatch::getMass(
  MASS_TYPE_T mass_type /// MONO or AVERAGE?
  ) {

  if (!mass_calculated_[mass_type]) {
    mass_[mass_type] = calcMass(mass_type);
    mass_calculated_[mass_type] = true;
  }
  return mass_[mass_type];
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

/**
 * Print one field in the tab-delimited output file, based on column index.
 */
void XLinkMatch::printOneMatchField(
  int      column_idx,             ///< Index of the column to print. -in
  MatchCollection* collection,  ///< collection holding this match -in 
  MatchFileWriter*    output_file,            ///< output stream -out
  int      scan_num,               ///< starting scan number -in
  FLOAT_T  spectrum_precursor_mz,  ///< m/z of spectrum precursor -in
  int      num_target_matches,            ///< target matches for this spectrum -in
  int      num_decoy_matches, ///< decoy matches for this spectrum -in
  int      b_y_total,              ///< total b/y ions -in
  int      b_y_matched             ///< Number of b/y ions matched. -in
) {

  carp(CARP_DETAILED_DEBUG, "XLinkMatch::printOneMatchField:%s", get_column_header(column_idx));

  switch ((MATCH_COLUMNS_T)column_idx) {

  case PEPTIDE_MASS_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, getMass(MONO));
    break;
  case PVALUE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, pvalue_);
    break;
  case SEQUENCE_COL:
    output_file->setColumnCurrentRow((MATCH_COLUMNS_T)column_idx, 
      getSequenceString());
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
    output_file->setColumnCurrentRow(
    (MATCH_COLUMNS_T)column_idx,
    getScore(XCORR_FIRST));
    break;
  case XCORR_SECOND_COL:
    output_file->setColumnCurrentRow(
    (MATCH_COLUMNS_T)column_idx,
    getScore(XCORR_SECOND));
    break;
  case PROTEIN_ID_X_COL:
    output_file->setColumnCurrentRow(
      (MATCH_COLUMNS_T)column_idx,
      getProteinIdXString());
    break;

  default:
    Match::printOneMatchField(column_idx,
      collection,
      output_file,
      scan_num,
      spectrum_precursor_mz,
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

