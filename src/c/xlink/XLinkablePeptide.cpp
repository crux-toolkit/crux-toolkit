/*************************************************************************//**
 * \file XLinkablePeptide.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for finding and defining the link sites on a peptide.
 ****************************************************************************/
#include "XLinkablePeptide.h"

#include "objects.h"
#include "modifications.h"

#include <iostream>

#include "XLink.h"

using namespace std;
using namespace Crux;

/**
 * Initialize object
 */
void XLinkablePeptide::init() {
  peptide_ = NULL;
  sequence_ = NULL;
  is_decoy_ = false;
  link_sites_.clear();
  
}

/**
 * Default constructor
 */
XLinkablePeptide::XLinkablePeptide() {
  init();
}

/**
 *  Constructor that defines the peptide as a sequence string
 */
XLinkablePeptide::XLinkablePeptide(
  char* sequence ///< the peptide sequence
  ) {
  init();
  sequence_ = sequence;
  peptide_ = NULL;
  is_decoy_ = false;
}

XLinkablePeptide::XLinkablePeptide(
  const XLinkablePeptide& xlinkablepeptide
  ) {
  init();
  peptide_ = xlinkablepeptide.peptide_->copyPtr();
  is_decoy_ = xlinkablepeptide.is_decoy_;
  link_sites_ = xlinkablepeptide.link_sites_;
}

XLinkablePeptide::XLinkablePeptide(
  XLinkablePeptide& xlinkablepeptide
) {
  init();
  peptide_ = xlinkablepeptide.peptide_->copyPtr();
  is_decoy_ = xlinkablepeptide.is_decoy_;
  link_sites_ = xlinkablepeptide.link_sites_;
}

/**
 * Constructor that defines the peptide and the linking sites
 */
XLinkablePeptide::XLinkablePeptide(
  Peptide* peptide, ///< the peptide object 
  vector<int>& link_sites ///< the linking sites
  ) {
  init();
  peptide_ = peptide->copyPtr();

  for (unsigned int idx=0;idx<link_sites.size();idx++) {
    addLinkSite(link_sites[idx]);
  }
  is_decoy_ = false;
}

/**
 * Constructor that generates the possible linking sites
 * on the peptide using the bondmap object
 */
XLinkablePeptide::XLinkablePeptide(
  Peptide* peptide, ///< the peptide object 
  XLinkBondMap& bondmap ///< the bond map
  ) {
  init();
  peptide_=peptide->copyPtr();
  findLinkSites(peptide_, bondmap, link_sites_);
}


/**
 * \returns whether a link at this index in the sequence
 * would prevent cleavage
 */
bool XLinkablePeptide::linkSeqPreventsCleavage(
  int seq_idx) {

  return linkSeqPreventsCleavage(peptide_, seq_idx);
}

bool XLinkablePeptide::linkSeqPreventsCleavage(
  Peptide* peptide,
  int seq_idx
  ) {

  if (seq_idx >= 0 && seq_idx < peptide->getLength()) {
 
    if (isCleavageSite(peptide, seq_idx)) {

      char aa = peptide->getSequencePointer()[seq_idx];

      const char* xlink_prevents_cleavage=get_string_parameter_pointer("xlink-prevents-cleavage");
      int idx = 0;
      while (xlink_prevents_cleavage[idx] != '\0') {
        if (xlink_prevents_cleavage[idx] == '*' || xlink_prevents_cleavage[idx] == aa) {
          return true;
        }
        idx++;
      }
    }
  }
  return false;
}

/**
 * \returns whether the link site prevents enzymatic
 * cleavage
 */
bool XLinkablePeptide::linkSitePreventsCleavage(
  size_t link_site_idx
  ) {

  if (link_site_idx >= numLinkSites()) {
    return false;
  }

  return linkSeqPreventsCleavage(getLinkSite(link_site_idx));
}


bool XLinkablePeptide::isMissedCleavageSite(
  Peptide* peptide,
  int seq_idx) {

  if (seq_idx == peptide->getLength()-1) {
    return false;
  }
  return isCleavageSite(peptide, seq_idx);
}
  
/**
 * \returns whether the peptide is cleavable at this index
 */
bool XLinkablePeptide::isCleavageSite(
  Peptide* peptide, ///< the peptide
  int seq_idx ///< the sequence index
  ) {

  //TODO - make this work for every enzyme
  if (seq_idx > peptide->getLength()-1) {
    return false;
  }

  char* sequence = peptide->getSequencePointer();

  if (sequence[seq_idx] == 'K' || sequence[seq_idx] == 'R') {
    if (seq_idx == peptide->getLength()-1) {
      return true;
    } else {
      return sequence[seq_idx+1] != 'P';
    }
  } else {
    return false;
  }
  
}

/**
 * given a peptide and a XLinkBondMap object,
 * generate the possible link sites
 */
void XLinkablePeptide::findLinkSites(
  Peptide* peptide,  ///< the peptide object -in
  XLinkBondMap& bondmap,  ///< the bond map -in 
  vector<int>& link_sites ///< the found link sites -out
  ) {

  int missed_cleavages = getMissedCleavageSites(peptide);
  int max_missed_cleavages = get_int_parameter("missed-cleavages")+
    2; // +2 because a self loop can prevent two cleavages from happening

  link_sites.clear();

  //find all modifications that can prevent xlink.
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  vector<AA_MOD_T*> prevents_xlink;

  for (int mod_idx = 0;mod_idx < total_mods; mod_idx++) {
    if (aa_mod_get_prevents_xlink(mod_list[mod_idx])) {
      prevents_xlink.push_back(mod_list[mod_idx]);
   }
  }

  //scan through the sequence and find linkable sites.
  MODIFIED_AA_T* mod_seq = peptide->getModifiedAASequence();
  int length = peptide->getLength();

  //scan the sequence for linkage sites.
  for (int seq_idx=0;seq_idx < length;seq_idx++) {
    if (bondmap.canLink(peptide, seq_idx)) {

      bool link_prevented = false;
      // 1st test, make sure that we don't link at the c-terminus
      // if the xlink prevents cleavage and the c-terminus is
      // a cleavage point.  
      if ((seq_idx == length-1) &&
           linkSeqPreventsCleavage(peptide, seq_idx)) {
        link_prevented = true;
      }

      //2nd test, make sure that there isn't a modification that
      //prevents the cross-link
      if (!link_prevented) {
        //check if the modification prevents xlink.
        for (unsigned int mod_idx=0;mod_idx<prevents_xlink.size();mod_idx++) {
          //if aa is modified by any modification that can prevent
          //cross-linking, then don't add the site as a link_site.
          if (is_aa_modified(mod_seq[seq_idx], prevents_xlink[mod_idx])) {
            link_prevented = true;
            break;
          }
        }
      }

      //3rd test make sure that we are within the allowed 
      //number of missed cleavages
      if (!link_prevented) {
        //if a link occurs here and prevents a cleavage, then don't
        //count the number of missed cleavages against the total number 
        //of missed cleavages
        int total_missed_cleavages = missed_cleavages;
        if (linkSeqPreventsCleavage(peptide, seq_idx)) {
          total_missed_cleavages--;
        }

        if (total_missed_cleavages > max_missed_cleavages) {
          link_prevented = true;
        }
      }
      //passes all three tests, this is a linkable site.
      if (!link_prevented) {
        //if it is a linkable site, then add it to the list.
        link_sites.push_back(seq_idx);
      }
    }
  }
  free(mod_seq);
}

/**
 * \returns the number of link sites on this peptide
 */
size_t XLinkablePeptide::numLinkSites() {
  return link_sites_.size();
}

/**
 * \returns whether this peptide is linkable or not
 */
bool XLinkablePeptide::isLinkable() {
  return numLinkSites() > 0;
}

void XLinkablePeptide::setDecoy(bool is_decoy) {

  is_decoy_ = is_decoy;
}

/**
 * \returns whether the peptide is a decoy or not
 */
bool XLinkablePeptide::isDecoy() {
  //return peptide_ -> isDecoy();
  return is_decoy_;
}

/**
 * \returns whether the peptide is linkable using the bond map
 */
bool XLinkablePeptide::isLinkable(
    Peptide* peptide, ///< the peptide object
    XLinkBondMap& bondmap ///< the bond map
  ) {

  vector<int> link_sites;
  findLinkSites(peptide, bondmap, link_sites);
  return link_sites.size() > 0;


}

/**
 * \returns the number of missed cleavages in the linkable peptide
 */
int XLinkablePeptide::getMissedCleavageSites() {
  return getMissedCleavageSites(peptide_);
}

int XLinkablePeptide::getMissedCleavageSites(
  Peptide* peptide ///< the peptide
) {
  //TODO - change this to reflect modifications that 
  //prevent cleavage (i.e. deadlinks).
  return peptide -> getMissedCleavageSites();
}

/**
 * \returns the sequence index of a link site
 */
int XLinkablePeptide::getLinkSite(
  int link_site_idx ///< the index of the link site
  ) {

  return link_sites_.at(link_site_idx);
}

/**
 * Adds a link site by sequence index
 * \returns the index of the added link site
 */
int XLinkablePeptide::addLinkSite(
  int seq_idx ///< the sequence index of the link site
  ) {
  link_sites_.push_back(seq_idx);
  return link_sites_.size()-1;
}

/**
 * \returns the peptide object associated with this XLinkablePeptide
 */
Peptide* XLinkablePeptide::getPeptide() {
  return peptide_;
}

/**
 * \returns the mass of the xlinkable peptide
 */
FLOAT_T XLinkablePeptide::getMass(
  MASS_TYPE_T mass_type
  ) const {
  if (peptide_ == NULL) {
    return Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
    return peptide_->calcModifiedMass(mass_type);
  }
}

/**
 * \returns an allocated sequence c-string.  Must be freed
 */
char* XLinkablePeptide::getSequence() {
  if (peptide_ == NULL) {
    return my_copy_string(sequence_);
  } else {
    char* seq = peptide_->getSequence();
    return seq;
  }
}

/**
 * \returns the allocated modified sequence. Must be freed.
 */
MODIFIED_AA_T* XLinkablePeptide::getModifiedSequence() {
  MODIFIED_AA_T* mod_seq = NULL;

  if (peptide_ == NULL) {
    convert_to_mod_aa_seq(sequence_, &mod_seq);
  } else {
    mod_seq = peptide_->getModifiedAASequence();
  }

  return mod_seq;
}

int findLink(vector<int>& link_sites, int link_site) {
  
  for (unsigned int idx=0;idx<link_sites.size();idx++) {
    if (link_sites[idx] == link_site)
      return idx;
  }
  return -1;

}

void doSwitch(vector<int>& link_sites, int siteA, int siteB) {
  int siteA_link_idx  = findLink(link_sites, siteA);
  int siteB_link_idx = findLink(link_sites, siteB);
  
  if (siteA_link_idx != -1 && siteB_link_idx != -1) {
    int temp_idx = link_sites[siteA_link_idx];
    link_sites[siteA_link_idx] = link_sites[siteB_link_idx];
    link_sites[siteB_link_idx] = temp_idx;
  } else if (siteA_link_idx != -1) {
    link_sites[siteA_link_idx] = siteB;
  } else if (siteB_link_idx != -1) {
    link_sites[siteB_link_idx] = siteA;
  }
}

static const int MAX_SHUFFLES = 5;

char* generateShuffledSequence(
  Peptide* peptide,
  vector<int>& link_sites
  ) {

  char* sequence = peptide->getSequence();

  int length = peptide->getLength();

  int num_shuffles = 0;
  do {
    
    //Don't move the n-term and c-term amino acids.
    int start_idx = 1;
    int end_idx = length - 2;

    while (start_idx < end_idx) {
      int switch_idx = get_random_number_interval(start_idx, end_idx);
      char temp_char = sequence[start_idx];
      sequence[start_idx] = sequence[switch_idx];
      sequence[switch_idx] = temp_char;
     
      doSwitch(link_sites, start_idx, switch_idx);

      ++start_idx;
    }


  } while (false /*equal_peptides(sequence, peptide)*/ && (num_shuffles < MAX_SHUFFLES));

  return sequence;


}

MODIFIED_AA_T* generateShuffledModSequence(
  Peptide* peptide,
  vector<int>& link_sites
  ) {

  MODIFIED_AA_T* sequence = peptide->getModifiedAASequence();
  int length = peptide->getLength();
  int start_idx = 0;
  int end_idx = length-1;
  int switch_idx = 0;
  MODIFIED_AA_T temp_aa = 0;

  // Do not move the first and last residue, regardless of enzyme
  ++start_idx;
  --end_idx;
  
  // shuffle from left to right, using the Knuth algorithm for shuffling.
  while(start_idx < end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_aa = sequence[start_idx];
    sequence[start_idx] = sequence[switch_idx];
    sequence[switch_idx] = temp_aa;
    doSwitch(link_sites, start_idx, switch_idx);
    ++start_idx;
  }

  return sequence;

}


XLinkablePeptide XLinkablePeptide::shuffle() {
  //cerr <<"XLinkablePeptide::shuffle():start"<<endl;
  Peptide* peptide = new Peptide(peptide_);
  XLink::addAllocatedPeptide(peptide);
  //cerr <<"Linksites"<<endl;
  vector<int> link_sites;
  for (unsigned int idx=0;idx<numLinkSites();idx++) {
    link_sites.push_back(getLinkSite(idx));
  }

  MODIFIED_AA_T* new_mod_seq = NULL;

  if(peptide->isModified()) {
    //cerr<<"Calling generateShuffledModSequence"<<endl;
    new_mod_seq = generateShuffledModSequence(peptide, link_sites);
    peptide->setDecoyModifiedSeq(new_mod_seq);
  } else {
    //cerr<<"Calling generateShuffledSequence"<<endl;
    char* new_seq =
      generateShuffledSequence(peptide, link_sites);
    convert_to_mod_aa_seq(new_seq, &new_mod_seq);
    peptide->setDecoyModifiedSeq(new_mod_seq);
    free(new_seq);
  }
  //cerr<<"Creating new linkable peptide"<<endl;
  XLinkablePeptide ans(peptide, link_sites);
  //cerr <<"XLinkablePeptide::shufle():done."<<endl;
  return ans;
} 

/*
string getAASequence() {
  char* seq = get_peptide_sequence(peptide_);
  string ans(seq);
  free(seq);
  return ans;
}
*/

/**
 * \returns the modified sequence string of the xlinkable peptide
 */
string XLinkablePeptide::getModifiedSequenceString() {
  if (peptide_ == NULL) {
    return string(sequence_);
  } else {
    char* seq = peptide_->getModifiedSequenceWithMasses(MOD_MASS_ONLY);
    string string_seq(seq);
    free(seq);
    return string_seq;
  }

}

bool XLinkablePeptide::isModified() {
  return peptide_->isModified();
}


bool compareXLinkablePeptideMass(
  const XLinkablePeptide& xpep1,
  const XLinkablePeptide& xpep2
) {

  return xpep1.getMass() < xpep2.getMass();
}
/**
 * \returns whether the peptide is less than (by lexical modified sequence)
 */
bool XLinkablePeptide::operator < (
    XLinkablePeptide other ///< the other XLinkablePeptide to compare to.
  ) const {
  
  XLinkablePeptide xlinkable = *this;
  return xlinkable.getModifiedSequenceString() < other.getModifiedSequenceString();

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
