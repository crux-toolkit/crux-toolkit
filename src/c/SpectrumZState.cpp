
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
}

SpectrumZState::SpectrumZState(
  FLOAT_T neutral_mass,
  int charge
  ) {

  neutral_mass_ = neutral_mass;
  charge_ = charge;
}

/** 
 * Default destructor
 */
SpectrumZState::~SpectrumZState() {

}

/**
 * \returns The charge for this z-state
 */
int SpectrumZState::getCharge() {
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
FLOAT_T SpectrumZState:: getSinglyChargedMass() {
  
  return neutral_mass_ + MASS_PROTON;
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
FLOAT_T SpectrumZState::getNeutralMass() {

  return neutral_mass_;
}

/**
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
