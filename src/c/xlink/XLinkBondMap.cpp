/*************************************************************************//**
 * \file XLinkBondMap.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for representing the potential cross-links for peptides.
 ****************************************************************************/

#include "XLinkBondMap.h"
#include "DelimitedFile.h"

using namespace std;

/**
 * Default constructor.
 */
XLinkBondMap::XLinkBondMap() {
}

/**
 * Constructor that initializes XLinkBondMap using the links string.
 * Format: A:B,C:D,... which means that a link can occur between residue
 * A and B, or C and D.
 */
XLinkBondMap::XLinkBondMap(
  string& links_string ///<link string
  ) {
  
  vector<string> bond_strings;

  DelimitedFile::tokenize(links_string, bond_strings, ',');

  for (unsigned int bond_idx = 0; bond_idx < bond_strings.size(); bond_idx++) {
    vector<string> link_site_strings;

    DelimitedFile::tokenize(bond_strings[bond_idx], link_site_strings, ':');

    if (link_site_strings.size() == 2) {
      XLinkSite site1(link_site_strings[0]);
      XLinkSite site2(link_site_strings[1]);
      (*this)[site1].insert(site2);
      (*this)[site2].insert(site1);
    } else {
      carp(CARP_FATAL,
        "bad format in %s when parsing %s",
        links_string.c_str(),
        bond_strings[bond_idx].c_str());
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

  for (XLinkBondMap::iterator iter1 = begin();
    iter1 != end(); ++iter1) {

    if (iter1->first.hasSite(peptide, idx1)) {
        for (set<XLinkSite>::iterator iter2 = iter1->second.begin();
          iter2 != iter1->second.end();
          ++iter2) {
      
        if (iter2->hasSite(peptide, idx2)) {
          return true;
        }
      }
    }
  }
  return false;
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
