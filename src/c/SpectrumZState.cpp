
/*************************************************************************//**
 * \file SpectrumZState.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE:  January 28, 2011
 * \brief code to support spectrum precursor z-states.
 ****************************************************************************/
#include "SpectrumZState.h"

#include "carp.h"
#include "mass.h"


/**
 * Default constructor
 */
SpectrumZState::SpectrumZState() {
  charge_ = 0;
  neutral_mass_ = 0;
  rtime_ = 0;
  area_ = 0;
}

SpectrumZState::SpectrumZState(
  FLOAT_T neutral_mass,
  int charge
  ) {

  neutral_mass_ = neutral_mass;
  charge_ = charge;
  rtime_ = 0;
  area_ = 0;
}

/**
 * copy constructor
 */
SpectrumZState::SpectrumZState(
  const SpectrumZState& other
) {
  charge_ = other.charge_;
  neutral_mass_ = other.neutral_mass_;
  rtime_ = other.rtime_;
  area_ = other.area_;
}

/** 
 * Default destructor
 */
SpectrumZState::~SpectrumZState() {

}

/**
 * \returns The charge for this z-state
 */
int SpectrumZState::getCharge() const {
  return charge_;
}

/**
 * Sets the neutral mass and charge for this z-state
 * using the m/z and charge.
 */

void SpectrumZState::setMZ(
  FLOAT_T mz,
  int charge
) {

  neutral_mass_ = (mz - MASS_PROTON) * (double)charge;
  charge_ = charge;

}

FLOAT_T SpectrumZState::getMZ() const {

  return (neutral_mass_ > 0) ?
    (neutral_mass_ + (FLOAT_T)charge_*MASS_PROTON) / (FLOAT_T)charge_ :
    0;
}


/**
 * Sets the neutral mass and charge for this z-state
 * using the singly charged mass.
 */
void SpectrumZState::setSinglyChargedMass(
  FLOAT_T mph,
  int charge
  ) {

  neutral_mass_ = mph - MASS_PROTON;
  charge_ = charge;
}

/**
 * \returns the m+h charged mass for this z-state
 */
FLOAT_T SpectrumZState:: getSinglyChargedMass() const {
  
  return (neutral_mass_ > 0) ?
    neutral_mass_ + MASS_PROTON :
    0;
}

/**
 * Sets the neutral mass for this z-state
 */
void SpectrumZState::setNeutralMass(
  FLOAT_T neutral_mass,
  int charge
  ) {

  neutral_mass_ = neutral_mass;
  charge_ = charge;
}

/**
 * \returns The neutral mass for this z-state
 */
FLOAT_T SpectrumZState::getNeutralMass() const {

  return (neutral_mass_ > 0) ? neutral_mass_ : 0;
}

/** 
 * Sets the retention time for this z-state
 */
void SpectrumZState::setRTime(
  FLOAT_T rtime
  ) {

  rtime_ = rtime;
}

/**
 * \returns The retention time for this z-state
 */
FLOAT_T SpectrumZState::getRTime() const {
  return rtime_;
}

/**
 * Sets the area for this z-state
 */
void SpectrumZState::setArea(
    FLOAT_T area
  ) {

  area_ = area;
}

/**
 * \returns The area for this z-state
 */
FLOAT_T SpectrumZState::getArea() const {
  return area_;
}


/**
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
