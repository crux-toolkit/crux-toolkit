/**
 * \file SpectrumZState.h
 * $Revision: 1.00 $
 * \brief Object for representing a MS2 spectrum's charge, and neutral mass.
 *******************************************************************************/
#ifndef SPECTRUMZSTATE_H
#define SPECTRUMZSTATE_H

#include "model/objects.h"
#include "util/utils.h"


/**
 * \class SpectrumZState
 * \brief Object for representing a spectrum precursor's 
 *  charge and neutral mass.
 */

class SpectrumZState {
 protected:
  int charge_;
  double neutral_mass_;

  /* EZ State fields */
  double rtime_;
  double area_;
 public:

  /**
   * Default constructor
   */
  SpectrumZState();

  /**
   * sets the neutral mass and charge
   */
  SpectrumZState(
    double neutral_mass, 
    int charge
  );

  /**
   * copy constructor
   */
  SpectrumZState(
    const SpectrumZState& other
  );

  /** 
   * Default destructor
   */
  virtual ~SpectrumZState();

  /**
   * \returns The charge for this z-state
   */
  int getCharge() const;

  
  /**
  * sets the m/z, charge for this z-state
  */
  void setMZ(
    double mz,
    int charge);

  double getMZ() const;


  /**
   * Sets the m+h charged mass for this z-state
   */
  void setSinglyChargedMass(
    double mph,
    int charge
    );

  /**
   * \returns the m+h charged mass for this z-state
   */
  double getSinglyChargedMass() const;

  /**
   * Sets the neutral mass for this z-state
   */
  void setNeutralMass(
    double neutral_mass,
    int charge
    );

  /**
   * \returns The neutral mass for this z-state
   */
  double getNeutralMass() const;
  
  /** 
   * Sets the retention time for this z-state
   */
  void setRTime(
    double rtime
  );

  /**
   * \returns The retention time for this z-state
   */
  double getRTime() const;

  /**
   * Sets the area for this z-state
   */
  void setArea(
    double area
  );

  /**
   * \returns The area for this z-state
   */
  double getArea() const;




};

/**
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


#endif

