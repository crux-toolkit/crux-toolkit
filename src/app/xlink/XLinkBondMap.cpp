/*************************************************************************
 * \file XLinkBondMap.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for representing the potential cross-links for peptides.
 ****************************************************************************/

#include "XLinkBondMap.h"
#include "util/Params.h"
#include "util/StringUtils.h"

using namespace std;
using namespace Crux;

/**
 * Default constructor.
 */
XLinkBondMap::XLinkBondMap() {
  ;
}

/**
 * Constructor that initializes XLinkBondMap using the links string.
 * Format: A:B,C:D,... which means that a link can occur between residue
 * A and B, or C and D.
 */
XLinkBondMap::XLinkBondMap(
  string& links_string ///<links string
  ) {
  init(links_string);
}

/**
 * Initializes XLinkBondMap using the links string.
 * Format: A:B,C:D,... which means that a link can occur between residue
 * A and B, or C and D.
 */
/**
 * link sites - Specify the two sets of amino acids that the cross-linker can 
 * connect. These are specified as two comma-separated sets of amino acids, 
 * with the two sets separated by a colon. Cross-links involving the terminus 
 * of a protein can be specified by using "nterm" or "cterm". For example, 
 * "K,nterm:Q" means that the cross linker can attach K to Q or the protein 
 * N-terminus to Q. Note that the vast majority of cross-linkers will 
 * operate on the following reactive groups: amine (K,nterm), 
 * carboxyl (D,E,cterm), sulfhydrl (C), acyl (Q) or amine+ (K,S,T,Y,nterm).
 */

void XLinkBondMap::init(
  string& links_string ///< links string
  ) {

  vector<string> bond_strings = StringUtils::Split(links_string, ':');
  if (bond_strings.size() != 2) {
    carp(CARP_FATAL, "invalid links string format:%s",links_string.c_str());
  }
  
  vector<string> link1_sites = StringUtils::Split(bond_strings[0], ',');
  vector<string> link2_sites = StringUtils::Split(bond_strings[1], ',');

  for (vector<string>::const_iterator idx1 = link1_sites.begin();
        idx1 != link1_sites.end();
       ++idx1) {
    XLinkSite site1(*idx1);
    for (vector<string>::const_iterator idx2 = link2_sites.begin();
	 idx2 != link2_sites.end();
	 ++idx2) {
      XLinkSite site2(*idx2);
      (*this)[site1].insert(site2);
      (*this)[site2].insert(site1);
    }
  }
}

/**
 * Default destructor
 */
XLinkBondMap::~XLinkBondMap() {
}

/**
 * \returns whether a cross-link can occur at a single position in the 
 * peptide (for deadlinks).
 */
bool XLinkBondMap::canLink(
  Peptide* peptide, ///<peptide object pointer
  int idx             ///<sequence index
  ) {
  for (XLinkBondMap::iterator iter = begin();
    iter != end(); ++iter) {

    if (iter->first.hasSite(peptide, idx)) {
      return true;
    }
  }
  return false;
}

/**
 * \returns whether a cross-link can occur between two positions in the 
 * peptide (for selfloops).
 */
bool XLinkBondMap::canLink(
    Peptide* peptide, ///<peptide object pointer
    int idx1,           ///<1st sequence idx
    int idx2            ///<2nd sequence idx
    ) {

  return canLink(peptide, peptide, idx1, idx2);
}

bool XLinkBondMap::canLink(
  XLinkablePeptide& pep,
  int link1_site,
  int link2_site
  ) {

  return canLink(
    pep.getPeptide(),
    pep.getLinkSite(link1_site),
    pep.getLinkSite(link2_site));
}

/**
 * \returns whether a cross-link can occur between two peptides at their 
 * respective sequence positions (for inter/intra links).
 */
bool XLinkBondMap::canLink(
  Peptide* peptide1,  ///<1st peptide object pointer 
  Peptide* peptide2,  ///<2nd peptide object pointer
  int idx1,             ///<1st peptide sequence idx
  int idx2              ///<2nd peptide sequence idx
  ) { //for inter/intra links

  for (XLinkBondMap::iterator iter1 = begin();
    iter1 != end(); ++iter1) {

    if (iter1->first.hasSite(peptide1, idx1)) {
      for (set<XLinkSite>::iterator iter2 = iter1->second.begin();
        iter2 != iter1->second.end();
        ++iter2) {
      
        if (iter2->hasSite(peptide2, idx2)) {
          return true;
        }
      }
    }
  }
  return false;
}

bool XLinkBondMap::canLink(
  XLinkablePeptide& pep1,
  XLinkablePeptide& pep2,
  int link1_site,
  int link2_site
  ) {
  
  return canLink(
    pep1.getPeptide(),
    pep2.getPeptide(),
    pep1.getLinkSite(link1_site),
    pep2.getLinkSite(link2_site));

}


bool XLinkBondMap::canLink(
  string& protein,
  int idx) {

  
  for (XLinkBondMap::iterator iter = begin();
    iter != end(); ++iter) {

    if (iter->first.hasSite(protein, idx)) {
      return true;
    }
  }
  return false;

}

/*                                                                                                                                                                                                                          
 * Local Variables:                                                                                                                                                                                                         
 * mode: c                                                                                                                                                                                                                  
 * c-basic-offset: 2                                                                                                                                                                                                        
 * End:                                                                                                                                                                                                                     
 */

