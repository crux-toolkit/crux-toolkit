/**
 * \file SelfLoopPeptide.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 18 September December 2014
 * \brief Object for Defining a self-loop peptide in an xlink search
 *****************************************************************************/
#include "SelfLoopPeptide.h"
#include "XLinkablePeptide.h"
#include "XLinkablePeptideIterator.h"
#include "XLinkPeptide.h"
#include "XLink.h"
#include "XLinkDatabase.h"
#include "util/GlobalParams.h"

#include "model/IonSeries.h"
#include "model/Ion.h"

#include <iostream>
#include <sstream>

using namespace std;

/**
 * Default Constructor
 */
SelfLoopPeptide::SelfLoopPeptide() {
  is_decoy_ = false;
}


/**
 * Constructor that defines a linkable peptide and the positions
 */
SelfLoopPeptide::SelfLoopPeptide(
  XLinkablePeptide& peptide, ///< linkable peptide
  int posA, ///< 1st link position on peptide
  int posB  ///< 2nd link position on peptide
  ) {
  
  is_decoy_ = false;
  linked_peptide_ = XLinkablePeptide(peptide);
  linked_peptide_.clearSites();
  linked_peptide_.addLinkSite(posA);
  linked_peptide_.addLinkSite(posB);

  link_pos_idx_.push_back(0);
  link_pos_idx_.push_back(1);

}

/**
 * Constructor that defines a linkable peptide and the positions
 */
SelfLoopPeptide::SelfLoopPeptide(
  char* peptide, ///< linkable peptide
  int posA, ///< 1st link position on peptide
  int posB ///< 2nd link position on peptide
  ) {

  linked_peptide_ = peptide;
  link_pos_idx_.push_back(posA);
  link_pos_idx_.push_back(posB);
  
}

/**
 * Adds Self-loop candidates to the collection
 */
void SelfLoopPeptide::addCandidates(
  FLOAT_T min_mass, ///< min mass
  FLOAT_T max_mass, ///< max mass
  bool is_decoy, ///< is decoy.
  XLinkMatchCollection& candidates ///< collection to add candidates
  ) {


  vector<SelfLoopPeptide>::iterator biter = XLinkDatabase::getSelfLoopBegin(is_decoy, min_mass);
  if (biter == XLinkDatabase::getSelfLoopEnd(is_decoy) ||
      biter -> getMass(GlobalParams::getIsotopicMass()) > max_mass) {
    return;
  } else {
    vector<SelfLoopPeptide>::iterator eiter = XLinkDatabase::getSelfLoopEnd(is_decoy);

    while (biter != eiter && biter->getMass(GlobalParams::getIsotopicMass()) <= max_mass) {
      biter->incrementPointerCount();
      candidates.add(&(*biter));
      ++biter;
    }
  }
}

/**
 * \returns the link position
 */
int SelfLoopPeptide::getLinkPos(
  int link_idx ///< link index (0 or 1)
  ) {
  return linked_peptide_.getLinkSite(link_pos_idx_[link_idx]);
}

/**
 * \returns self-loop candidate
 */
XLINKMATCH_TYPE_T SelfLoopPeptide::getCandidateType() {
  return SELFLOOP_CANDIDATE;
}

/**
 * \returns the sequence string, marking the link sites with (X,Y).
 */
string SelfLoopPeptide::getSequenceString() {

  string seq = linked_peptide_.getModifiedSequenceString();
  ostringstream oss;
  oss << seq << " (" << (getLinkPos(0)+1) << "," << (getLinkPos(1)+1) << ")";
  string svalue = oss.str();

  return svalue;
}

/**
 * \returns the mass of the self-loop peptide
 */
FLOAT_T SelfLoopPeptide::calcMass(
  MASS_TYPE_T mass_type ///< AVERAGE or MONO
  ) {
  
  return linked_peptide_.getMass(mass_type) + XLinkPeptide::getLinkerMass();
}

/**
 * \returns a shuffled version of self-loop candidate
 */
void SelfLoopPeptide::shuffle(vector<XLinkMatch*>& decoys) {
  carp(CARP_DEBUG, "SelfLoopPeptide::shuffle");
  SelfLoopPeptide* decoy = new SelfLoopPeptide();

  decoy->linked_peptide_ = linked_peptide_.shuffle();
  decoy->link_pos_idx_.push_back(link_pos_idx_[0]);
  decoy->link_pos_idx_.push_back(link_pos_idx_[1]);
  decoy->is_decoy_ = true;
  decoy->target_ = this;
  decoys.push_back(decoy);
}

SelfLoopPeptide* SelfLoopPeptide::getUnshuffledTarget() {
  return(target_);  
}

std::string SelfLoopPeptide::getUnshuffledSequence() {
  
  if (is_decoy_) {
    return(getUnshuffledTarget()->getSequenceString());
  } else {
    return(getSequenceString());
  }
}

/**
 * Predicts the ions for the self loop candidate
 */
void SelfLoopPeptide::predictIons(
  IonSeries* ion_series, ///< ion vector to place ions in
  int charge ///< charge state
  ) {
  
  const char* seq = linked_peptide_.getSequence();
  const MODIFIED_AA_T* mod_seq = linked_peptide_.getModifiedSequencePtr();
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();

  unsigned int first_site = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));
  unsigned int N = strlen(seq);


  //iterate through the ions and modify the ones that have the linker 
  //attached.
  vector<Ion*> to_remove;

  for (IonIterator ion_iter = ion_series->begin();
    ion_iter != ion_series->end();
    ++ion_iter) {

    Ion* ion = *ion_iter;

    bool keep_ion = false;
    bool modify_ion = false;
    unsigned int cleavage_idx = ion->getCleavageIdx();
   
    if (ion->isForwardType()) {
      if (cleavage_idx > first_site) {
        if (cleavage_idx <= second_site) {
          keep_ion = false;
        } else {
          keep_ion = true;
          modify_ion = true;
        }
      } else {
        keep_ion = true;
        modify_ion = false;
      }
    } else {
      if (cleavage_idx >= (N-second_site)) {
        if (cleavage_idx < (N-first_site)) {
          keep_ion = false;
        } else {
          keep_ion = true;
          modify_ion = true;
        }
      } else {
        keep_ion = true;
        modify_ion = false;
      }
    }

    if (keep_ion) {
      if (modify_ion) {
        FLOAT_T mass_z = ion->getMassZ();
        int charge = ion->getCharge();
        double mass = (mass_z -MASS_PROTON) * (double)charge;
        mass += XLinkPeptide::getLinkerMass();
        mass_z = (mass + MASS_PROTON * (double)charge) / (double)charge;
        ion->setMassZ(mass_z);
        if (isnan(ion->getMassZ())) {
          carp(CARP_FATAL, "NAN5");
        }
        
      }
    } else {
      to_remove.push_back(ion);
    }
  }

  for (unsigned int idx=0;idx < to_remove.size();idx++) {
    ion_series->removeIon(to_remove[idx]);
    delete to_remove[idx];
  }

}

/**
 * \returns the sequence of ion
 */
string SelfLoopPeptide::getIonSequence(
  Ion* ion ///< ion
  ) {

  string ion_sequence = ion->getPeptideSequence();

  unsigned int cleavage_idx = ion->getCleavageIdx();

  unsigned int first_site  = min(getLinkPos(0), getLinkPos(1));
  unsigned int second_site = max(getLinkPos(0), getLinkPos(1));

  bool is_linked = false;
  if (ion->isForwardType()) {
    is_linked = (cleavage_idx > first_site);
  } else {
    is_linked = (cleavage_idx >= (ion_sequence.length() - second_site));
  }

  string subseq;

  //cerr<<"creating substring"<<endl;
  if (ion->isForwardType()) {
    subseq = ion_sequence.substr(0, cleavage_idx);
  } else {
    subseq = ion_sequence.substr(ion_sequence.length() - cleavage_idx, ion_sequence.length());
  }

  if (is_linked) {
    return subseq + string("*");
  } else {
    return subseq;
  }
}

/**
 *\returns the peptide object of self-loop peptide (0 is only valid)
 */
Crux::Peptide* SelfLoopPeptide::getPeptide(
  int peptide_idx
  ) {
  if (peptide_idx == 0) {
    return linked_peptide_.getPeptide();
  } else {
    return NULL;
  }
}

/**
 *\returns the number of missed cleavages
 */
int SelfLoopPeptide::getNumMissedCleavages() {
  char missed_cleavage_link_site = 'K';

  int link1_site = getLinkPos(0);
  int link2_site = getLinkPos(1);

  set<int> skip;

  Crux::Peptide* pep = linked_peptide_.getPeptide();

  char* seq = pep->getSequencePointer();

  if (seq[link1_site] == missed_cleavage_link_site) {
    skip.insert(link1_site);
  }

  if (seq[link2_site] == missed_cleavage_link_site) {
    skip.insert(link2_site);
  }

  return pep->getMissedCleavageSites(skip);

}

/**
 *\returns whether the peptide is modified by a variable mod
 */
bool SelfLoopPeptide::isModified() {

  return linked_peptide_.isModified();
}

/**
 *\returns "target" or "decoy"
 */
string SelfLoopPeptide::getDecoyType() {
  string returnValue = "target";
  if (is_decoy_) {
    returnValue = "decoy";
  }
  return(returnValue);
}

bool compareSelfLoopPeptideMass(
				const SelfLoopPeptide& spep1,
				const SelfLoopPeptide& spep2) {
  return spep1.getMassConst(GlobalParams::getIsotopicMass()) < spep2.getMassConst(GlobalParams::getIsotopicMass());

}

bool compareSelfLoopPeptideMassToFLOAT(
				       const SelfLoopPeptide& spep1,
				       FLOAT_T mass) {

  return spep1.getMassConst(GlobalParams::getIsotopicMass()) < mass;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
