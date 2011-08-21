/*************************************************************************//**
 * \file XLinkSite.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for representing cross-linkable site on a peptide.
 ****************************************************************************/

#include "XLinkSite.h"

using namespace std;

/**
 * Default constructor
 */
XLinkSite::XLinkSite() {
  type_ = XLINKSITE_UNKNOWN;
  aa_ = 0;
}

/**
 * Constructor that determines the type from the supplied site string 
 */
XLinkSite::XLinkSite(
  string& site_string ///<string describing the cross-linkable site
  ) {

  if (site_string == "nterm") {

    type_ = XLINKSITE_NTERM;

  } else if (site_string == "*") {

    type_ = XLINKSITE_ALL;

  } else if (site_string.length() == 1) {

    type_ = XLINKSITE_AA;
    aa_=site_string[0];

  } else {

    carp(CARP_FATAL, "Unknown link site:%s",site_string.c_str());
  }
}

/**
 * Default destructor
 */
XLinkSite::~XLinkSite() {
}

/**
 * \returns whether the peptide contains this site at the supplied sequence index.
 */
bool XLinkSite::hasSite(
  Peptide* peptide, ///<peptide object pointer 
  int idx             ///<sequence index
  ) const {
  
  switch (type_) {
    case XLINKSITE_NTERM:
      if (idx == 0) {
        
        PeptideSrc* src = peptide->getPeptideSrc();
        while (src != NULL) {
          carp(CARP_DEBUG,"nterm peptide start_idx:%d",src->getStartIdx());
          if (src->getStartIdx() == 1) {
            return true;
          }
          src = src->getNextAssociation();
        }
      }
      return false;
      break;
    case XLINKSITE_ALL:
      return true;
      break;
    case XLINKSITE_AA:
      return peptide->getSequencePointer()[idx] == aa_;
      break;
    case XLINKSITE_UNKNOWN:
    default:
      carp(CARP_FATAL, "Xlink site not set!");
  };
  return false;
}

/**
 * \returns whether this xlinksite is equal to the passed in xlinksite
 */
bool XLinkSite::operator == (
  XLinkSite& xlink_site_obj ///<XLinkSite to compare to
  ) {

  XLinkSite site1 = *this;
  XLinkSite site2 = xlink_site_obj;

  if (site1.type_ == site2.type_) {
    if (site1.type_ == XLINKSITE_AA) {
      return site1.aa_ == site2.aa_;
    }
    return true;
  } else {
    return false;
  }
}

/**
 * \returns whether this xlinksite is less than the passed in xlinksite
 */
bool XLinkSite::operator < (
    const XLinkSite& xlink_site_obj ///<XLinkSite to compare to
  ) const {

  XLinkSite site1 = *this;
  XLinkSite site2 = xlink_site_obj;

  if (site1.type_ == site2.type_) {
    if (site1.type_ == XLINKSITE_AA) {
      return site1.aa_ < site2.aa_;
    }
    return false;
  } else {
    return site1.type_ < site2.type_;
  }
}

