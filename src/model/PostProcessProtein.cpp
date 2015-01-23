/**
 * \file PostProcessProtein.cpp
 * AUTHOR : Sean McIlwain
 * $Revision: 1.00 $
 * \brief Object for representing a protein that doesn't exist within the database
 */

#include "PostProcessProtein.h"
#include <iostream>

using namespace std;

  
/**
 * Default constructor
 */
PostProcessProtein::PostProcessProtein() : Protein() {
}

/**
 * Default destructor
 */
PostProcessProtein::~PostProcessProtein() {
}


/**
  * \returns the index of the sequence in the sequences_ vector and
  * returns the index.  If sequence is not found, then add to the
  * list of sequences and return the last index
  */
int PostProcessProtein::findStart(
  string sequence, ///< the sequence to find
  string prev_aa,  ///< the previous aa of the sequence
  string next_aa   ///< the next aa of the sequence
  ) {

  int ans = -1;

  for (int idx = 0; idx < sequences_.size();idx++) {

    if (sequences_[idx] == sequence) {
      ans = idx+1;
      break;
    }
  }

  if (ans == -1) {
    sequences_.push_back(sequence);
    prev_aas_.push_back(prev_aa);
    next_aas_.push_back(next_aa);
    ans = sequences_.size();
  }

  return ans;

}

/**
 * \returns the ith sequence
 */
char* PostProcessProtein::getSequence(
  int offset ///< The offset (or sequence index) for the sequence
  ) {

  char* ans = my_copy_string(sequences_[offset].c_str());

  return ans;
}

/**
 * \returns the ith sequence pointer
 */
char* PostProcessProtein::getSequencePointer(
  int offset ///< The offset (or sequence index) for the sequence
  ) {
  return (char*)(sequences_[offset].c_str());

}

char PostProcessProtein::getNTermFlankingAA(
  int offset ///< The offset (or sequence index) for the flanking AA
  ) {
  char ans = prev_aas_.at(offset)[0];
  if (ans == '\0') {
    carp_once(CARP_WARNING, "Missing nterm flanking for protein:%s offset:%d",
      getIdPointer(), offset); 
    //just return - for now.
      ans = '-';
  }
  carp(CARP_DETAILED_DEBUG, "protein:%s offset:%d pflank:%c(%d)", 
    getIdPointer(), offset, ans,(int)ans);

  return ans;
}

char PostProcessProtein::getCTermFlankingAA(
  int offset ///< The offset (or sequence index) for the flanking AA
  ) {

  char ans = next_aas_.at(offset)[0];
  if (ans == '\0') {
    carp_once(CARP_WARNING, "Missing cterm flanking aa for protein:%s offset:%d",getIdPointer(), offset);
    ans = '-';
  }
  carp(CARP_DETAILED_DEBUG, "protein:%s offset:%d nflank:%c(%d)", getIdPointer(), offset, ans,(int)ans);

  return ans;
}

/**
 * \returns true indicating that this is a PostProcessProtein object
 */
bool PostProcessProtein::isPostProcess() {

  return true;
}

/**
 * ends the program since we can't get the protein length from this
 * post process protein object
 */
unsigned int PostProcessProtein::getLength() {

  if (sequence_ == NULL) {
    carp_once(CARP_WARNING, "Need protein sequence in order to calculate protein length.\n"
                     "   Please provide protein fasta or index using the protein-database parameter\n"
                     "   Protein %s doesn't have the full sequence", getIdPointer());

    return 0;
  }

  return length_;

}

Database* PostProcessProtein::getDatabase() {
  return NULL;
}

