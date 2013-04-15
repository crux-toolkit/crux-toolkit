/**
 * \file LinkedPeptide.cpp 
 * AUTHOR: Sean McIlwain and Paul Draghicescu
 * CREATE DATE: 18 January 2012
 * \brief Object for keeping track of a (non-)crosslinked peptide.
 *****************************************************************************/
//TODO Change cout/cerrs to CARP.

#include "LinkedPeptide.h"
#include "XHHC_Peptide.h"

using namespace Crux;

/**
 * mass of the linker
 */
FLOAT_T LinkedPeptide::linker_mass_;

/**
 * Initializes the object
 */
void LinkedPeptide::init() {

  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  charge_ = 0;
  type_ = NUMBER_ION_TYPES;
  mass_[MONO] = 0.0;
  mass_[AVERAGE] = 0.0;
  mz_[MONO] = 0.0;
  mz_[AVERAGE] = 0.0;

}

/**
 * \returns a blank LinkedPeptide object
 */
LinkedPeptide::LinkedPeptide() {
  init();
}

/**
 * \returns a linked peptide object.
 * constructor for a linked peptide. If A or B null, then
 * a self link will be created. If an index is -1, a link to nothing
 * will be created.
 */
LinkedPeptide::LinkedPeptide(
  char* peptide_A, ///< First peptide
  char* peptide_B, ///< Second peptide (can be NULL)
  int posA, ///< Index of link on first peptide
  int posB, ///< Index of link on second peptide
  int charge ///< charge of product
  ) {

  init();
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
  charge_ = charge;
  decoy_ = false;
  type_ = (ION_TYPE_T)NULL;
  XHHC_Peptide pepA = XHHC_Peptide(peptide_A);
  // if a self link or dead end
  if (peptide_B == NULL) {
     pepA.addLink(posA);
    if (posB >= 0)
      pepA.addLink(posB);
    peptides_.push_back(pepA);
  } else {
    carp(CARP_DETAILED_DEBUG, "adding links at %d and %d", posA, posB);
    XHHC_Peptide pepB = XHHC_Peptide(peptide_B);
    pepA.addLink(posA);
    pepB.addLink(posB);
    peptides_.push_back(pepA);
    peptides_.push_back(pepB);
  }
  getMZ(MONO);
  getMZ(AVERAGE);
}

/**
 * \returns a blank LinkedPeptide object with 
 * the charge set
 */
LinkedPeptide::LinkedPeptide(
    int charge ///< charge of product
  ) {

  init();
  charge_ = charge;
  decoy_ = false;
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
}

/**
 * Destructor
 */
LinkedPeptide::~LinkedPeptide() {
}

/**
 * Sets the linker_mass_ static variable
 */
void LinkedPeptide::setLinkerMass(
  FLOAT_T linker_mass ///< the linker mass
  ) {

  linker_mass_ = linker_mass;
}

/**
 * \returns the linker_mass_ static variable
 */
FLOAT_T LinkedPeptide::getLinkerMass() {
  return linker_mass_;
}

/**
 * \returns a reference to the internal
 * peptides vector
 */
vector<XHHC_Peptide>& LinkedPeptide::getPeptides() {

  return peptides_;
}

/**
 * \returns the charge of the LinkedPeptide
 */
int LinkedPeptide::getCharge() {
  return charge_;
}

/**
 * sets the charge of the LinkedPeptide
 */
void LinkedPeptide::setCharge(
  int charge ///< the charge of the LinkedPeptide
  ) {

  charge_ = charge;
}

/**
 * sets the IonType for the LinkedPeptide
 */
void LinkedPeptide::setIonType(
  ION_TYPE_T type ///< The ion type
  ) {
  
  type_ = type;
}

/**
 * \returns the ion type for this LinkedPeptide
 */
ION_TYPE_T LinkedPeptide::getIonType() {
  return type_;
}

/**
 * \returns the number of peptides in this linked
 * Peptide (should be 1 or 2)
 */
int LinkedPeptide::size() {
  return peptides_.size();
}

/**
 * indicates that this LinkedPeptide is a
 * decoy
 */
void LinkedPeptide::setDecoy() {
  decoy_ = true;
}

/**
 * \returns whether the LinkedPeptide is a
 * decoy
 */
bool LinkedPeptide::isDecoy() {
  return decoy_;
}

/**
 * adds a XHHC_Peptide to the list of peptide
 * in this LinkedPeptide object
 */
void LinkedPeptide::addPeptide(
  XHHC_Peptide& peptide
  ) {

  peptides_.push_back(peptide);
  mass_calculated_[MONO] = false;
  mass_calculated_[AVERAGE] = false;
}

/**
 * \returns whether the LinkedPeptide is a
 * crosslinked peptide
 */
bool LinkedPeptide::isCrossLinked() {
  return size() == 2;
}

/**
 * \returns whether the Peptide is a normal/linear peptide
 */
bool LinkedPeptide::isLinear() {
  return (peptides_.size() == 1 && peptides_[0].linkSite() == -1);
} 

/**
 * \returns whether the Peptide is a deadend
 */
bool LinkedPeptide::isDeadEnd() {
  return (peptides_.size() == 1 && peptides_[0].getNumLinks() == 1);
}

/**
 * \returns whether the Peptide is a self-loop
 */
bool LinkedPeptide::isSelfLoop() {
  return (peptides_.size() == 1 && peptides_[0].getNumLinks() == 2);
}

/**
 * /returns the mass of the LinkedPeptide
 */
FLOAT_T LinkedPeptide::getMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {

  if (!mass_calculated_[mass_type])
    calculateMass(mass_type);
  return mass_[mass_type];
   
}

/**
 * calculates the mass of the LinkedPeptide
 */
void LinkedPeptide::calculateMass(
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ){

  // remove H2O from mass if it's a b-ion

  mass_[mass_type] = Peptide::calcSequenceMass((char*)peptides_[0].getSequence().c_str(), mass_type);   

  if (peptides_[0].getNumLinks() > 0) {
    mass_[mass_type] += linker_mass_;
  }

  if (size() == 2) {
    mass_[mass_type] += Peptide::calcSequenceMass((char*)peptides_[1].getSequence().c_str(), mass_type);
  }

  if (type_ == B_ION) {
    if (mass_type == MONO) {
      mass_[mass_type] -= MASS_H2O_MONO;
    }
    else {
      mass_[mass_type] -= MASS_H2O_AVERAGE;
    }
  } 
  mass_calculated_[mass_type] = true;
}

/**
 * \returns the m/z of the LinkedPeptide
 */
FLOAT_T LinkedPeptide::getMZ (
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  )  {

  if (mass_type == MONO) {
    mz_[MONO] = (getMass(MONO) + MASS_PROTON * charge_) / charge_;
  } else {
    mz_[AVERAGE] = (getMass(AVERAGE) + MASS_H_AVERAGE * charge_) / charge_;
  }
  return mz_[mass_type];
}

/**
 * splits the LinkedPeptide into b-y ions
 */
void LinkedPeptide::split(
  vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///< the b-y ion pairs
  ) {

  // split between every amino acid on every peptide in the
  // linked peptide.
  bool is_loop = false;
  XHHC_Peptide peptideA = peptides_[0];  
  XHHC_Peptide peptideB = peptides_[0];
  // dead end
  if (isDeadEnd()) {
   peptideB.setSequence("");
  }
  // if a loop
  if (peptideA.getNumLinks() == 2) {
    is_loop = true;
  } else if (isCrossLinked()) {
    peptideB = peptides_[1];
  }
  for (int idx = 1; idx < peptideA.getLength(); ++idx) {
    peptideA.splitAt(idx, ion_pairs, charge_, peptideB, is_loop);
  }
 
  if (isCrossLinked()) {
    for (int idx = 1; idx < peptideB.getLength(); ++idx) {
      peptideB.splitAt(idx, ion_pairs, charge_, peptideA, is_loop);
    }
  } 
}

/**
 * splits the first LinkedPeptide into b-y ions
 */
void LinkedPeptide::splitA(
  vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///< the b-y ion pair
  ) {

  XHHC_Peptide peptideA = peptides_[0];
  XHHC_Peptide peptideB = peptides_[1];

  for (int idx = 1; idx < peptideA.getLength(); ++idx) {
    peptideA.splitAt(idx, ion_pairs, charge_, peptideB, false);
  }
}

/**
 * splits the second LinkedPeptide into b-y ions
 */
void LinkedPeptide::splitB(
  vector<pair<LinkedPeptide, LinkedPeptide> >& ion_pairs ///< the b-y ion pair
  ) {

  XHHC_Peptide peptideA = peptides_[0];
  XHHC_Peptide peptideB = peptides_[1];

  for (int idx = 1; idx < peptideB.getLength(); ++idx) {
    peptideB.splitAt(idx, ion_pairs, charge_, peptideA, false);
  }
}

/**
 * used for sorting LinkedPeptides by mass
 */
bool operator < (
  const LinkedPeptide &lp1, 
  const LinkedPeptide &lp2) {

  return lp1.mass_ < lp2.mass_;
  //return lp1.mz < lp2.mz;
}

/**
 * prints the LinkedPeptide to a stream
 */
std::ostream &operator<< (std::ostream& os, LinkedPeptide& lp) {
  vector<XHHC_Peptide> peptides = lp.getPeptides();
  ostringstream link_positions;
  link_positions << "(";
  for (int idx = 0; idx < peptides[0].getLength(); ++idx) {
	if (peptides[0].hasLinkAt(idx)) link_positions << (idx+1) << "," ;
  }
  if (peptides.size() == 2) {
    for (int idx = 0; idx < peptides[1].getLength(); ++idx) {
	if (peptides[1].hasLinkAt(idx)) link_positions << (idx+1) << ")";
    }
    return os << peptides[0].getSequence() << "," << peptides[1].getSequence() << link_positions.str();// << " +" << lp.charge();
  }
  string final = link_positions.str();
  if (final.length() > 1) final.erase(final.end()-1);
  return os << peptides[0].getSequence() << final << ")";// +" << lp.charge();
}

bool compareMassAverage(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {
  if (lp1.mass_calculated_[AVERAGE] && lp2.mass_calculated_[AVERAGE]) { 
    return lp1.mass_[AVERAGE] < lp2.mass_[AVERAGE];
  } else {
    carp(CARP_FATAL, "LinkedPeptide Average mass not calculated!");
    return false;
  }
}

bool compareMassMono(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {
  if (lp1.mass_calculated_[MONO] && lp2.mass_calculated_[MONO]) {
    return lp1.mass_[MONO] < lp2.mass_[MONO];
  } else {
    carp(CARP_FATAL, "LinkedPeptide Mono mass not calculated!");
    return false;
  }
}

bool compareMZAverage(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {

  return lp1.mz_[AVERAGE] < lp2.mz_[AVERAGE];
}

bool compareMZMono(const LinkedPeptide& lp1, const LinkedPeptide& lp2) {
  return lp1.mz_[MONO] < lp2.mz_[MONO];
}


/**
 * Sorts a vector of LinkedPeptide by mass
 */
void LinkedPeptide::sortByMass(
  vector<LinkedPeptide>& linked_peptides, ///< the LinkedPeptides to sort 
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {
  //TODO : should we put code here to make sure that we have
  //calculated the all of the masses for a particular mass type?

  if (mass_type == MONO) {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMassMono);
  }
  else {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMassAverage);
  }
}

/**
 * Sorts a vector of LinkedPeptide by m/z
 */
/**
 * Sorts a vector of LinkedPeptide by mass
 */
void LinkedPeptide::sortByMZ(
  vector<LinkedPeptide>& linked_peptides, ///< the LinkedPeptides to sort 
  MASS_TYPE_T mass_type ///< MONO or AVERAGE
  ) {

  IF_CARP(CARP_DEBUG,
    //sanity check
    for (unsigned int idx=0;idx < linked_peptides.size();idx++) {
      if (!linked_peptides.at(idx).mass_calculated_[mass_type]) {
        carp(CARP_ERROR, "mass not calculated for LinkedPeptide!");
        linked_peptides.at(idx).getMZ(mass_type);
      }
    }
  )

  if (mass_type == MONO) {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMZMono);
  }
  else {
    sort(linked_peptides.begin(), linked_peptides.end(), compareMZAverage);
  }
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
