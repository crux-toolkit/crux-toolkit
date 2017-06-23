/*************************************************************************
 * \file XLinkablePeptide.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  Febuary 22, 2011
 * \brief  Object for finding and defining the link sites on a peptide.
 ****************************************************************************/
#include "XLinkablePeptide.h"

#include "XLinkIonSeriesCache.h"

#include "model/objects.h"
#include "util/modifications.h"
#include "util/Params.h"
#include "util/GlobalParams.h"
#include "model/Ion.h"
#include "model/IonSeries.h"

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
  xcorr_link_idx_ = -1;
  predict_ions_call_count_ = 0;
  index_ = -1;
  mod_seq_ = NULL;
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
  sequence_ = my_copy_string(sequence);
  peptide_ = NULL;
  is_decoy_ = false;
}

XLinkablePeptide::XLinkablePeptide(
  const XLinkablePeptide& xlinkablepeptide
  ) {
  init();
  if (xlinkablepeptide.peptide_) {
    peptide_ = xlinkablepeptide.peptide_;
  } else if (sequence_) {
    sequence_ = my_copy_string(xlinkablepeptide.sequence_);
  }
  is_decoy_ = xlinkablepeptide.is_decoy_;
  link_sites_ = xlinkablepeptide.link_sites_;
  xcorr_link_idx_ = xlinkablepeptide.xcorr_link_idx_;
  xcorr_ = xlinkablepeptide.xcorr_;
  index_ = xlinkablepeptide.index_;
  for (size_t idx=0;idx<NUMBER_MASS_TYPES;idx++) {
    if (xlinkablepeptide.mass_calculated_[idx]) {
      mass_[idx] = xlinkablepeptide.mass_[idx];
      mass_calculated_[idx] = true;
    } else {
      mass_calculated_[idx] = false;
    }
  }
}

XLinkablePeptide::XLinkablePeptide(
  XLinkablePeptide& xlinkablepeptide
) {
  init();
  peptide_ = xlinkablepeptide.peptide_;
  is_decoy_ = xlinkablepeptide.is_decoy_;
  link_sites_ = xlinkablepeptide.link_sites_;
  xcorr_link_idx_ = xlinkablepeptide.xcorr_link_idx_;
  xcorr_ = xlinkablepeptide.xcorr_;
  index_ = xlinkablepeptide.index_;
    for (size_t idx=0;idx<NUMBER_MASS_TYPES;idx++) {
    if (xlinkablepeptide.mass_calculated_[idx]) {
      mass_[idx] = xlinkablepeptide.mass_[idx];
      mass_calculated_[idx] = true;
    } else {
      mass_calculated_[idx] = false;
    }
  }  


}

/**
 * Constructor that defines the peptide and the linking sites
 */
XLinkablePeptide::XLinkablePeptide(
  Peptide* peptide, ///< the peptide object 
  vector<int>& link_sites ///< the linking sites
  ) {
  init();
  peptide_ = peptide;

  for (unsigned int idx = 0; idx < link_sites.size(); idx++) {
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
  XLinkBondMap& bondmap, ///< the bond map
  int additional_cleavages
  ) {
  init();
  peptide_ = peptide;
  findLinkSites(peptide_, bondmap, link_sites_, additional_cleavages);
}

XLinkablePeptide::~XLinkablePeptide() {
  //cerr <<"XLinkablePeptide::~XLinkablePeptide"<<endl;
  if (mod_seq_ != NULL) {
    freeModSeq(mod_seq_);
  }
  if (sequence_) {
    std::free(sequence_);
  }
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

      const string& xlink_prevents_cleavage = GlobalParams::getXLinkPreventsCleavage();
      for (string::const_iterator i = xlink_prevents_cleavage.begin();
           i != xlink_prevents_cleavage.end();
           i++) {
        if (*i == '*' || *i == aa) {
          return true;
        }
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
  vector<int>& link_sites, ///< the found link sites -out
  int additional_cleavages ///< 0 for xlink, 1 for self-loops
  ) {

  int missed_cleavages = getMissedCleavageSites(peptide);
  int max_missed_cleavages = GlobalParams::getMissedCleavages() +
    additional_cleavages; // +1 because a self loop can prevent two cleavages from happening

  link_sites.clear();

  //find all modifications that can prevent xlink.
  AA_MOD_T** mod_list = NULL;
  int total_mods = get_all_aa_mod_list(&mod_list);
  vector<AA_MOD_T*> prevents_xlink;

  for (int mod_idx = 0;mod_idx < total_mods; mod_idx++) {
    if (mod_list[mod_idx]->getPreventsXLink()) {
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
        for (unsigned int mod_idx = 0; mod_idx < prevents_xlink.size(); mod_idx++) {
          //if aa is modified by any modification that can prevent
          //cross-linking, then don't add the site as a link_site.
          if (prevents_xlink[mod_idx]->isModified(mod_seq[seq_idx])) {
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
  freeModSeq(mod_seq);
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


void XLinkablePeptide::clearSites() {
  link_sites_.clear();
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
FLOAT_T XLinkablePeptide::calcMass(
  MASS_TYPE_T mass_type
  ) {
  if (peptide_ == NULL) {
    return Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
    return peptide_->calcModifiedMass(mass_type);
  }
}

/**
 * \returns an allocated sequence c-string.  Must be freed
 */
const char* XLinkablePeptide::getSequence() {
  if (peptide_ == NULL) {
    return(sequence_);
  } else {
    return(peptide_->getSequence());
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
const MODIFIED_AA_T* XLinkablePeptide::getModifiedSequencePtr() {
  if (mod_seq_ == NULL) {
    mod_seq_ = getModifiedSequence();
  } else {
  }
  return(mod_seq_);
}

int findLink(vector<int>& link_sites, int link_site) {
  
  for (unsigned int idx = 0; idx < link_sites.size(); idx++) {
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
  while(start_idx < end_idx) {
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
  for (unsigned int idx = 0; idx < numLinkSites(); idx++) {
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
  for (size_t idx=0;idx<NUMBER_MASS_TYPES;idx++) {
    if (mass_calculated_[idx]) {
      ans.mass_calculated_[idx] = true;
      ans.mass_[idx] = mass_[idx];
    }
  }
  ans.is_decoy_ = true;
  return ans;
} 

XLinkablePeptide* XLinkablePeptide::getCachedDecoy() {
  if (decoy_ == NULL) {
    
    decoy_ = new XLinkablePeptide(shuffle());
  }
  return decoy_;
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
    return peptide_->getModifiedSequenceWithMasses();
  }
}

bool XLinkablePeptide::isModified() {
  return peptide_->isModified();
}

void XLinkablePeptide::setXCorr(size_t link_idx, FLOAT_T xcorr) {

  xcorr_link_idx_ = link_idx;
  xcorr_ = xcorr;

}

FLOAT_T XLinkablePeptide::getXCorr() const {
  if (xcorr_link_idx_ != -1) {
    return xcorr_;
  }
  carp(CARP_FATAL, "Xcorr not set!");
  return 0;

}

bool XLinkablePeptide::hasXCorr() const {
  //return (xcorr_link_idx_ != -1);
  return false;
}


bool compareXLinkableXCorr(
  const XLinkablePeptide& xpep1,
  const XLinkablePeptide& xpep2
  ) {
  return xpep1.getXCorr() > xpep2.getXCorr();
}

bool compareXLinkableXCorrPtr(
  const XLinkablePeptide* xpep1,
  const XLinkablePeptide* xpep2
  ) {
  return xpep1 -> getXCorr() > xpep2 -> getXCorr();
}



bool compareXLinkablePeptideMass(
  const XLinkablePeptide& xpep1,
  const XLinkablePeptide& xpep2
) {
  bool ret = xpep1.getMassConst(GlobalParams::getIsotopicMass()) < xpep2.getMassConst(GlobalParams::getIsotopicMass());
  return(ret);
}

bool compareXLinkablePeptideMassToFLOAT(
					const XLinkablePeptide& xpep1,
					const FLOAT_T& mass) {
  return xpep1.getMassConst(GlobalParams::getIsotopicMass()) < mass;
}

bool compareXLinkablePeptideMassToFLOAT2(
					 const FLOAT_T& mass,
					 const XLinkablePeptide& xpep1) {
  return xpep1.getMassConst(GlobalParams::getIsotopicMass()) > mass;
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

void XLinkablePeptide::predictIons(
  IonSeries* ion_series,
  int charge,
  int link_idx,
  FLOAT_T mod_mass,
  bool clear
  ) {
  const char* seq = NULL;
  IonSeries* cached_ions = NULL;
  bool cached = false;
  if (GlobalParams::getXLinkUseIonCache()) {
    cached_ions = XLinkIonSeriesCache::getXLinkablePeptideIonSeries(*this, charge);
    bool cached = cached_ions != NULL;
  }
  if (!cached) {
    cached_ions = new IonSeries(ion_series->getIonConstraint(), charge);
    seq = getSequence();
    cached_ions->update(seq, getModifiedSequencePtr());
    cached_ions->predictIons();
  }
  int link_pos = link_sites_.at(link_idx);
  int seq_len;
  if (peptide_) {
    seq_len = peptide_->getLength();
  } else {
    seq_len = strlen(sequence_);
  }
  if (clear) {
    ion_series->clear();
  }
  //modify the necessary ions and add to the ion_series
  IonIterator eiter = cached_ions->end();
  for (IonIterator ion_iter = cached_ions->begin(); 
    ion_iter != eiter; 
    ++ion_iter) { 
    Ion* src_ion = *ion_iter; 
    Ion* ion = src_ion;
    unsigned int cleavage_idx = ion->getCleavageIdx(); 
    if (ion->isForwardType()) { 
      if (cleavage_idx > (unsigned int)link_pos) {
        ion = Ion::newIon();
	if (seq) {
	  Ion::copy(src_ion, ion, seq);
	} else {
          Ion::copy(src_ion, ion, sequence_);
	}
        FLOAT_T mass = ion->getMassFromMassZ() + mod_mass;
        ion->setMassZFromMass(mass); 
        if (isnan(ion->getMassZ())) { 
          carp(CARP_FATAL, "NAN3"); 
        } 
      } 
    } else { 
      if (cleavage_idx >= (seq_len-(unsigned int)link_pos)) { 
        ion = Ion::newIon();
	if (seq) {
	  Ion::copy(src_ion, ion, seq);
	} else {
	  Ion::copy(src_ion, ion, sequence_);
	}
        FLOAT_T mass = ion->getMassFromMassZ() + mod_mass;
        ion->setMassZFromMass(mass); 
        if (isnan(ion->getMassZ())) { 
          carp(CARP_FATAL, "NAN4"); 
        } 
      } 
    }
    ion_series->addIon(ion);
  }
  
  if (!cached) {
    delete cached_ions;
  }
  if (seq) {
    std::free((char*)seq);
  }
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
