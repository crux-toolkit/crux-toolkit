/**
 * \file SpectrumZState.h
 * $Revision: 1.00 $
 * \brief Object for representing a MS2 spectrum's charge, and neutral mass.
 *******************************************************************************/
#ifndef SPECTRUMZSTATE_H
#define SPECTRUMZSTATE_H

#include "objects.h"
#include "utils.h"


/**
 * \class SpectrumZState
 * \brief Object for representing a spectrum precursor's 
 *  charge and neutral mass.
 */

class SpectrumZState {
 protected:
  int charge_;
  FLOAT_T neutral_mass_;

 public:

  /**
   * Default constructor
   */
  SpectrumZState();

  /**
   * sets the neutral mass and charge
   */
  SpectrumZState(
    FLOAT_T neutral_mass, 
    int charge
  );

  /** 
   * Default destructor
   */
  virtual ~SpectrumZState();

  /**
   * \returns The charge for this z-state
   */
  int getCharge();

  
  /**
  * sets the m/z, charge for this z-state
  */
  void setMZ(
    FLOAT_T mz,
    int charge);


  /**
   * Sets the m+h charged mass for this z-state
   */
  void setSinglyChargedMass(
    FLOAT_T mph,
    int charge
    );

  /**
   * \returns the m+h charged mass for this z-state
   */
  FLOAT_T getSinglyChargedMass();

  /**
   * Sets the neutral mass for this z-state
   */
  void setNeutralMass(
    FLOAT_T neutral_mass,
    int charge
    );

  /**
   * \returns The neutral mass for this z-state
   */
  FLOAT_T getNeutralMass();
  





};

/**
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


#endif

