/**
 * \file IonConstraint.h 
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * \brief Object for defining constraints on an ion series.
 *****************************************************************************/
#ifndef IONCONSTRAINT_H
#define IONCONSTRAINT_H

#include "objects.h"

#include "IonSeries.h"
#include "Ion.h"

/**
 * \class ion_constraint
 * \brief An object to represent the contraints which the ions in this
 * series obey.
 */
class IonConstraint {
 protected:
  bool use_neutral_losses_; ///< Should ions include neutral losses
  int modifications_[MAX_MODIFICATIONS]; 
    ///< an array to indicate which modifications to perform
  MASS_TYPE_T mass_type_; 
    ///< the mass_type to use MONO|AVERAGE
  int max_charge_; 
    ///< maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type_; 
    ///< the ion types the peptide series should include
  bool precursor_ion_; 
    ///< does a precursor-ion satisfy this constraint
  int min_charge_; 
  bool exact_modifications_; 
    ///< TRUE  = ints in modfications array indicate exact number of mods
    ///< FALSE = ints in modfications array indicate maximum number of mods
  unsigned int pointer_count_; ///< Number of pointers referencing me 

  /**
   *Initializes an IonConstraint object, called from constructor 
   */
  void init();

 public:

  /**
   *\returns an empty heap allocated ion_constraint
   */

  IonConstraint();

  /**
   * modification, all modifications 0
   * add more modifications as needed using the set_ion_constraint_modification
   *\returns a new heap allocated ion_constraint
   */
  IonConstraint(
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    ION_TYPE_T ion_type, ///< the ion types the peptide series should include
    bool precursor_ion  ///< should include precursor ion?
  );

  /**
   * \brief Create a new ion constraint based on the score type and the
   * charge of the peptide to be modeled.  Uses other
     * new_ion_constraint_ methods for some types.
   *
   * \returns A newly allocated ion constraint.
   */
  static IonConstraint* newIonConstraintSmart(
    SCORER_TYPE_T score_type,
    int charge
    );

  /**
   * modification, sets all fields for sequest settings
   *\returns a new heap allocated ion_constraint
   */
  static IonConstraint* newIonConstraintSequest(
    MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
    int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    ION_TYPE_T ion_type, ///< the ion types the peptide series should include
    bool precursor_ion ///< should include precursor ion?
    );

  /**
   * modification, sets all fields for GMTK settings
   *\returns a new heap allocated ion_constraint
   */
  static IonConstraint* newIonConstraintGmtk(
    int charge ///< the charge of the peptide for which to predict ions
  );

  /**
   * modification, sets all fields for sequest Sp scoring settings
   *\returns a new heap allocated ion_constraint
   */
  static IonConstraint* newIonConstraintSequestSp(
    int max_charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    );

  /**
   * modification, sets all fields for Sequest Xcorr scoring settings
   * make B, Y, A type ions
   *\returns a new heap allocated ion_constraint
   */
  static IonConstraint* newIonConstraintSequestXcorr(
    int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
    );

  /**
   * Frees an allocated ion_constraint object.
   */
  static void free(IonConstraint* ion_constraint);

  virtual ~IonConstraint();

  /**
   * copies the ion_constraint pointer
   */
  static IonConstraint* copy(
    IonConstraint* ion_constraint
  );

  /**
   * copies ion_constraint object from src to dest
   * must pass in a memory allocated ION_CONSTRAINT_T dest
   */
  static void copy(
    IonConstraint* src,///< ion_constraint to copy from -in
    IonConstraint* dest///< ion_constraint to copy to -out
    );

  /** 
   * Determines if a ion satisfies a ion_constraint.
   * \returns TRUE if the constraint is satisified. FALSE if not.
   */
  bool isSatisfied(
    Ion* ion ///< query ion -in
    );

  /**
   * \returns ION_TYPE for this constraint
   */
  ION_TYPE_T getIonType();


  /**
   * Sets the modification count
   * can only add isotopes
   */
  void setModification(
    ION_MODIFICATION_T mod_type, ///< ion modification type -in
    int count  ///< the count of the modification -in  
    );

  /**
   * Sets the exact modification boolean to exactness criteria
   * and if exactness is true sets min_charge = max_charge.
   * In other words, the constraint is now exact, in that it refers to a
   * particular ion series, charge states, and modification state, as opposed
   * to e.g. b-ions of charge state +1 or +2, or with or without NH3 loss
   */
  void setExactness(
    bool exactness ///< whether to be exact or not -in
    );

  /**
   * gets the modification count for specific mod_type
   */
  int getModification(
    ION_MODIFICATION_T mod_type ///< ion modification type -in
    );

  int* getModifications();

  /**
   * gets the mass type of the ion_constraint
   */
  MASS_TYPE_T getMassType();

  /**
   * get the maximum charge of the IonConstraint
   */
  int getMaxCharge();


  bool getPrecursorIon();

  /**
   * get the neutral loss field of the ion constraint.
   */
  bool getUseNeutralLosses();

};

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
#endif
