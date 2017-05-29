
#include "MonoLinkPeptide.h"
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
MonoLinkPeptide::MonoLinkPeptide() : LinearPeptide() {
}

/**
 * Constructor with sequence
 */
MonoLinkPeptide::MonoLinkPeptide(
  char* sequence ///< sequence string
  ) : LinearPeptide(sequence) {
}

/**
 * Constructor from a Crux Peptide
 */
MonoLinkPeptide::MonoLinkPeptide(
  Crux::Peptide* peptide ///< peptide object
  ) : LinearPeptide(peptide) {
}

MonoLinkPeptide::~MonoLinkPeptide() {
}


/**
 *Add candidates to the XLinkMatchCollection that are mono-links
 */
void MonoLinkPeptide::addCandidates(
  FLOAT_T min_mass,  ///< min mass
  FLOAT_T max_mass,  ///< max mass
  bool is_decoy, ///< generate decoy canidates
  XLinkMatchCollection& candidates ///< Vector of candidate -inout
  ) {

  vector<MonoLinkPeptide>::iterator siter = XLinkDatabase::getMonoLinkBegin(is_decoy, min_mass);
  if (siter == XLinkDatabase::getMonoLinkEnd(is_decoy) ||
      siter->getMass(GlobalParams::getIsotopicMass()) > max_mass) {
    return;
  } else {
    vector<MonoLinkPeptide>::iterator eiter = XLinkDatabase::getMonoLinkEnd(is_decoy, siter, max_mass);

    while (siter != eiter && siter->getMass() <= max_mass) {
      siter->incrementPointerCount();
      MonoLinkPeptide& lpeptide = *siter;
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
XLINKMATCH_TYPE_T MonoLinkPeptide::getCandidateType() {
  return DEADLINK_CANDIDATE;
}

/**
 * \returns the sequence of the peptide in string format
 */
string MonoLinkPeptide::getSequenceString() {
  ostringstream oss;
  if (peptide_ == NULL) {
    oss << sequence_;

  } else {
    oss << peptide_->getModifiedSequenceWithMasses();
  }

  //oss << " ()";
  return oss.str();

}


bool compareMonoLinkPeptideMass(
  const MonoLinkPeptide& pep1, 
  const MonoLinkPeptide& pep2) {

  return pep1.getMassConst(GlobalParams::getIsotopicMass()) < pep2.getMassConst(GlobalParams::getIsotopicMass());

}
bool compareMonoLinkPeptideMassToFLOAT(
  const MonoLinkPeptide& pep1,
  FLOAT_T mass
  ) {
  return pep1.getMassConst(GlobalParams::getIsotopicMass()) < mass;
}

bool compareFLOATToMonoLinkPeptideMass(
  const FLOAT_T& mass,
  const MonoLinkPeptide& pep1) {
  return pep1.getMassConst(GlobalParams::getIsotopicMass()) > mass;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
