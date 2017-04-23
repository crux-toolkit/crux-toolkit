
#include "LinearPeptide.h"
#include "XLink.h"
#include "model/ModifiedPeptidesIterator.h"
#include "model/IonSeries.h"
#include "util/GlobalParams.h"
#include "XLinkDatabase.h"

#include <iostream>

using namespace std;

/**
 * Default constructor
 */
LinearPeptide::LinearPeptide() : XLinkMatch() {
  peptide_ = NULL;
  sequence_ = NULL;
}

/**
 * Constructor with sequence
 */
LinearPeptide::LinearPeptide(
  char* sequence ///< sequence string
  ) : XLinkMatch() {
  peptide_ = NULL;
  sequence_ = sequence;
}

/**
 * Constructor from a Crux Peptide
 */
LinearPeptide::LinearPeptide(
  Crux::Peptide* peptide ///< peptide object
  ) : XLinkMatch() {
  peptide_ = peptide;
  sequence_ = NULL;
  this->setNullPeptide(peptide_->isDecoy());
}

LinearPeptide::~LinearPeptide() {
}


/**
 *Add candidates to the XLinkMatchCollection that are linear
 */
void LinearPeptide::addCandidates(
  FLOAT_T min_mass,  ///< min mass
  FLOAT_T max_mass,  ///< max mass
  bool is_decoy, ///< generate decoy canidates
  XLinkMatchCollection& candidates ///< Vector of candidate -inout
  ) {

  vector<LinearPeptide>::iterator siter = XLinkDatabase::getLinearBegin(is_decoy, min_mass);
  if (siter == XLinkDatabase::getLinearEnd(is_decoy) ||
      siter->getMass(GlobalParams::getIsotopicMass()) > max_mass) {
    return;
  } else {
    vector<LinearPeptide>::iterator eiter = XLinkDatabase::getLinearEnd(is_decoy, siter, max_mass);

    while (siter != eiter && siter->getMass() <= max_mass) {
      siter->incrementPointerCount();
      LinearPeptide& lpeptide = *siter;
      if (lpeptide.getMass() < min_mass || lpeptide.getMass() > max_mass) {
        carp(CARP_DEBUG,
             "The mass %g of peptide %s is outside the precursor range of %g-%g.",
             lpeptide.getMass(), lpeptide.getSequenceString().c_str(),
             min_mass, max_mass);
        return;
      } else {
        //carp(CARP_INFO, "Add linear candidate");
        candidates.add(&(*siter));
        ++siter;
      }
    }
  }
}

/**
 * returns the candidate type, either a deadlink or a linear candidate
 */
XLINKMATCH_TYPE_T LinearPeptide::getCandidateType() {
  return LINEAR_CANDIDATE;
}

/**
 * \returns the sequence of the peptide in string format
 */
string LinearPeptide::getSequenceString() {
  ostringstream oss;
  if (peptide_ == NULL) {
    oss << sequence_;

  } else {
    oss << peptide_->getModifiedSequenceWithMasses();
  }

  //oss << " ()";
  return oss.str();

}

/**
 * \returns the mass of the peptide
 */
FLOAT_T LinearPeptide::calcMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {
  
  if (peptide_ == NULL) {
    return Crux::Peptide::calcSequenceMass(sequence_, mass_type);
  } else {
    return peptide_->calcModifiedMass(mass_type);
  }
}

/**
 *\returns a shuffled version of the peptide
 */
void LinearPeptide::shuffle(vector<XLinkMatch*>& decoys) {
  string seq = getSequenceString();
  Crux::Peptide* decoy_peptide = new Crux::Peptide(peptide_);

  decoy_peptide->transformToDecoy();

  XLink::addAllocatedPeptide(decoy_peptide);

  LinearPeptide* decoy = new LinearPeptide(decoy_peptide);
  
  string decoy_seq = decoy->getSequenceString();
  //  carp(CARP_INFO, "Shuffled: %s", decoy_seq.c_str());
  
  decoy->setNullPeptide(true);
  decoys.push_back(decoy);
}

/**
 * predicts the ions for this peptide
 */
void LinearPeptide::predictIons(
  IonSeries* ion_series, ///< ion series to place the ions
  int charge ///< charge state of the peptide
  ) {

  const char* seq = NULL;
  MODIFIED_AA_T* mod_seq = NULL;
  if (peptide_ == NULL) {
    seq = sequence_;
    convert_to_mod_aa_seq(seq, &mod_seq); 
  } else {
    seq = peptide_->getSequence();
    mod_seq = peptide_->getModifiedAASequence();
  }
  
  //carp(CARP_INFO, "predicting ions for :%s", seq);
  /*
  char* mseq = peptide_->getModifiedSequenceWithMasses(MOD_MASSES_SEPARATE);
  carp(CARP_INFO, "mod seq:%s", mseq);
  free(mseq);
  */
  ion_series->setCharge(charge);
  ion_series->update(seq, mod_seq);
  ion_series->predictIons();

  freeModSeq(mod_seq);

}

/**
 *\returns the ion sequence as a string
 */
string LinearPeptide::getIonSequence(
  Ion* ion ///< ion object
  ) {

  string seq_str = string(sequence_);

  int cleavage_idx = ion->getCleavageIdx();
  if (ion->isForwardType() == B_ION) {
    return seq_str.substr(0, cleavage_idx);
  } else {
    return seq_str.substr(seq_str.length()-cleavage_idx, seq_str.length());
  }
}

/**
 * \returns the peptide for this match
 */
Crux::Peptide* LinearPeptide::getPeptide(int peptide_idx) {
  if (peptide_idx == 0) {
    return peptide_;
  } else {
    return NULL;
  }
}

/**
 *\returns the number of missed cleavages
 */
int LinearPeptide::getNumMissedCleavages() {
  set<int> skip;
  return peptide_->getMissedCleavageSites(skip);
}

/**
 * \returns whether this peptide is modified by a variable mod
 */
bool LinearPeptide::isModified() {

  return peptide_->isModified();
}

/**
 *\returns "target" or "decoy"
 */
std::string LinearPeptide::getDecoyType() {
  string returnValue = "target";
  if (is_decoy_) {
    returnValue = "decoy";
  }
  return(returnValue);
}

bool compareLinearPeptideMass(
  const LinearPeptide& pep1, 
  const LinearPeptide& pep2) {

  return pep1.getMassConst(GlobalParams::getIsotopicMass()) < pep2.getMassConst(GlobalParams::getIsotopicMass());

}
bool compareLinearPeptideMassToFLOAT(
  const LinearPeptide& pep1,
  FLOAT_T mass
  ) {
  return pep1.getMassConst(GlobalParams::getIsotopicMass()) < mass;
}

bool compareFLOATToLinearPeptideMass(
  const FLOAT_T& mass,
  const LinearPeptide& pep1) {
  return pep1.getMassConst(GlobalParams::getIsotopicMass()) > mass;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
