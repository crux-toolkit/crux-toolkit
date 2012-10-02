/**
 * \file XHHC_Peptide.cpp 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 19 January 2012
 * \brief Object for keeping track of the links on a
 * cross-linkable peptide.
 *****************************************************************************/
#include "XHHC_Peptide.h"
#include "LinkedPeptide.h"
#include "crux-utils.h"


using namespace std;
using namespace Crux;
/**
 * \returns a blank XHHC_Peptide object
 */
XHHC_Peptide::XHHC_Peptide() {
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
}

/**
 * \returns XHHC_Peptide object with the sequence initialized
 */
XHHC_Peptide::XHHC_Peptide(
  string sequence
  ) {

  num_links_ = 0;
  sequence_ = sequence;
  length_ = sequence.length();
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  for (int idx = 0;idx < length_;idx++) {
    links_.push_back(false);
  }
}

/**
 * Destructor
 */
XHHC_Peptide::~XHHC_Peptide() {
}

/**
 * \returns whether a link exists at the index
 */
bool XHHC_Peptide::hasLinkAt(
  int index ///< amino acid index
  ) {
  
  return links_[index];
}

/**
 * \returns the sequence length of the peptide
 */
int XHHC_Peptide::getLength() {
  return length_;
}

/**
 * returns whether a link exists
 */
bool XHHC_Peptide::hasLink() {
  return num_links_ > 0;
}

/**
 * \returns the peptide sequence
 */
string XHHC_Peptide::getSequence() {
  return sequence_;
}

/**
 * \returns whether we have initialized the peptide
 */
bool XHHC_Peptide::isEmpty() {
  return length_ == 0;
}

/**
 * \returns the number of links
 */
int XHHC_Peptide::getNumLinks() {
  return num_links_;
}

/**
 * sets the peptide sequence
 */
void XHHC_Peptide::setSequence(
  string sequence ///< the sequence of the peptide
  ) {
  sequence_ = sequence;
  length_ = sequence.length();  
}

/**
 * removes a link from the peptide
 */
void XHHC_Peptide::removeLink(
  int index ///< the amino acid index
  ) {

  links_[index] = false;
  num_links_--;
}

/**
 * \returns the first index of a link, -1 if none exist
 */
int XHHC_Peptide::linkSite() {

  for (int idx = 0; idx < length_; ++idx) {
    if (hasLinkAt(idx)) 
      return idx;
  }
  return -1;
}

/**
 * \returns the mass of the peptide
 */
FLOAT_T XHHC_Peptide::getMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {

  if (mass_calculated_[mass_type]) 
    return mass_[mass_type];
  else {
    mass_[mass_type] = Peptide::calcSequenceMass((char*)sequence_.c_str(), mass_type);
    mass_calculated_[mass_type] = true;
    return mass_[mass_type];
  }
}

/**
 * add a link at a amino acid index
 */
void XHHC_Peptide::addLink(
  int index ///< the amino acid index
  ) {
  links_[index] = true;
  num_links_++;
}

/**
 * cleaves the peptide at an index, adding b and y ions
 * skips if a cleave site on self loop between the link
 */
void XHHC_Peptide::splitAt(
  int index, ///< the index to cleave
  vector<pair<LinkedPeptide, LinkedPeptide> >& pairs, ///< the b-y ion pairs
  int charge, ///< charge of the peptide
  XHHC_Peptide& other, ///< the other peptide we are linked to
  bool is_loop ///< Is this a self loop?
  ) {

  bool self_flag = false;
  // for every charge state
  for (int c = 0; c <= charge; ++c) {
    XHHC_Peptide pepA = XHHC_Peptide(sequence_.substr(0, index));
    XHHC_Peptide pepB = XHHC_Peptide(sequence_.substr(index, length_ - index));
    LinkedPeptide linkedA = LinkedPeptide(c);
    LinkedPeptide linkedB = LinkedPeptide(charge - c);
    self_flag = true;
    // for every position on pepA
    for (int idx = 0; idx < index; idx++) {
      if (hasLinkAt(idx)) {
        pepA.addLink(idx); 
        if (!other.isEmpty() && !is_loop) linkedA.addPeptide(other);
        // if a loop, skip cleavages in between link sites (same mass as precursor)
	if (is_loop) self_flag = !self_flag;
      }
    }
    if (!self_flag) continue;
    // for every position on pepB
    for (int idx = index; idx < length_; idx++) {
      if (hasLinkAt(idx)) {
        pepB.addLink(idx - index);
	if (!other.isEmpty() && !is_loop) linkedB.addPeptide(other);
	//else (self_flag = !self_flag);
      }
    } 
    linkedA.addPeptide(pepA);
    linkedB.addPeptide(pepB);
    pairs.push_back(pair<LinkedPeptide, LinkedPeptide> (linkedA, linkedB));
  }
}

/**
 * \returns a shuffled peptide, preserving any links
 */
XHHC_Peptide XHHC_Peptide::shuffle() {
  string shuffled = string(sequence_);
  XHHC_Peptide shuffled_peptide = XHHC_Peptide(sequence_);
  for (size_t idx = 0; idx < shuffled.length(); ++idx) {
    if (hasLinkAt(idx)) shuffled_peptide.addLink(idx);
  }

  int start_idx = 1;
  int end_idx = length_ - 2;
  int switch_idx = 0;
  char temp_char = 0;
  while(start_idx <= end_idx){
    switch_idx = get_random_number_interval(start_idx, end_idx);
    temp_char = shuffled[start_idx];
    shuffled[start_idx] = shuffled[switch_idx];
    shuffled[switch_idx] = temp_char;
    if (shuffled_peptide.hasLinkAt(switch_idx)) {
      //if not a self loop
      if (!shuffled_peptide.hasLinkAt(start_idx)) {
        shuffled_peptide.removeLink(switch_idx);
        shuffled_peptide.addLink(start_idx);
      }
    } else if (shuffled_peptide.hasLinkAt(start_idx)) {
      shuffled_peptide.removeLink(start_idx);
      shuffled_peptide.addLink(switch_idx);
    }
    ++start_idx;
  }
  shuffled_peptide.setSequence(shuffled);
  return shuffled_peptide;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
