/*************************************************************************//**
 * \file XLinkSite.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for representing cross-linkable site on a peptide.
 ****************************************************************************/

#include "XLinkSite.h"

using namespace std;
using namespace Crux;

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

  if (site_string == "cterm") {

    type_ = XLINKSITE_CTERM;

  } else if (site_string == "nterm") {

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
    case XLINKSITE_CTERM:
      
      if (idx == peptide->getLength() - 1) {
        //cerr <<"Cterm peptide:"<<peptide->getSequence()<<" "<<idx<<":"<<(int)peptide->getLength()<<endl;
        vector<PeptideSrc*>& srcs = peptide->getPeptideSrcVector();
        for (size_t idx2 = 0 ; idx2 < srcs.size() ; idx2++ ) {
          int start_idx = srcs[idx2]->getStartIdx();
          //cerr << "start:"<<start_idx<<" "<<(start_idx+idx)<<" "<<srcs[idx2]->getParentProtein()->getLength()<<endl;


          if ((idx + start_idx) >= srcs[idx2]->getParentProtein()->getLength()) {
            return true;
          }
        }
      }
      return false;
      break;
    case XLINKSITE_NTERM:
      if (idx == 0) {
        vector<PeptideSrc*>& srcs = peptide->getPeptideSrcVector();
        for (size_t idx = 0; idx < srcs.size(); idx++ ) {
          if (srcs[idx]->getStartIdx() == 1) {
            return true;
          }
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

bool XLinkSite::hasSite(
  string& protein_sequence,
  int idx) const {

  switch(type_) {
    case XLINKSITE_CTERM:
      return idx == protein_sequence.length() - 1;
    case XLINKSITE_NTERM:
      return idx == 0;
      break;
    case XLINKSITE_ALL:
      return true;
      break;
    case XLINKSITE_AA:
      return protein_sequence.at(idx) == aa_;
      break;
    case XLINKSITE_UNKNOWN:
    default:
      carp(CARP_FATAL, "XLink site not set!");
  }
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

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */
