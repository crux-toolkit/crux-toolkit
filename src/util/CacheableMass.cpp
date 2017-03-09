/*************************************************************************//**
 * \file CacheableMass.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 28 Oct 2016
 * \brief abstract object to Store the calculated mass of an instantied object
 ****************************************************************************/

#include "CacheableMass.h"
#include "io/carp.h"

/**
 * Initialized the object
 */
void CacheableMass::init() {
  for (size_t idx=0;idx < NUMBER_MASS_TYPES;idx++) {
    mass_calculated_[idx] = false;
  }
}

/**
 * Constructor
 */
CacheableMass::CacheableMass() {
  init();
}

/**
 * Destructor
 */
CacheableMass::~CacheableMass() {
}

/*
 * If the mass is already calculated for this mass type, return it
 * Otherwise call calcMass to set it and return the result.
 */
FLOAT_T CacheableMass::getMass(
  MASS_TYPE_T mass_type ///< Mass type
  ) {
 
  if (!mass_calculated_[mass_type]) {
    mass_[mass_type] = calcMass(mass_type);
    mass_calculated_[mass_type] = true;
  }

  return mass_[mass_type];

}

FLOAT_T CacheableMass::getMassConst(
  MASS_TYPE_T mass_type
  ) const {

  if (!mass_calculated_[mass_type]) {
    carp(CARP_FATAL, "Internal Error: mass not calculated yet! %d", mass_type);
  }
  return (mass_[mass_type]);
}

/*
 * Copies the calculated mass arrray
 */
void CacheableMass::copy(
  CacheableMass* src,
  CacheableMass* dest
  ) {

  for (size_t idx=0;idx < NUMBER_MASS_TYPES;idx++) {
    if (src->mass_calculated_[idx]) {
      dest->mass_calculated_[idx] = true;
      dest->mass_[idx] = src->mass_[idx];
    } else {
      dest->mass_calculated_[idx] = false;
    }
  }

}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
