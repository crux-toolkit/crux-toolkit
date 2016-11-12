/**
 * \file CacheableMass.h 
 * \brief abstract object to Store the calculated mass of an instantied object
 */

/*
 * AUTHOR: Sean McIlwain
 * CREATE DATE: Oct. 28, 2016
 * $Revision: 1.00 $
 *****************************************************************************/
#ifndef CACHEABLEMASS_H_
#define CACHEABLEMASS_H_

#include "model/objects.h"

class CacheableMass {

 protected:
  bool mass_calculated_[NUMBER_MASS_TYPES]; ///< is mass calculated?
  FLOAT_T mass_[NUMBER_MASS_TYPES]; ///< calculated mass

 public:
  /*
   * Initializes the object
   */
  void init();
  
  CacheableMass();
  
  virtual ~CacheableMass();

  /*
   * Override this method to return the calculated mass of the object
   */
  virtual FLOAT_T calcMass(
    MASS_TYPE_T mass_type
    )=0;


  /*
   * If the mass is already calculated for this mass type, return it
   * Otherwise call calcMass to set it and return the result.
   */
  FLOAT_T getMass(
    MASS_TYPE_T mass_type = MONO ///< Type of mass to return
    );

  FLOAT_T getMassConst(
    MASS_TYPE_T mass_type = MONO ///< Type of mass to return
  ) const;

  /*
   * Copies the calculated mass arrray
   */
  static void copy(
    CacheableMass* src, ///< Source 
	CacheableMass* dest ///< Dest
	);
 
};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
