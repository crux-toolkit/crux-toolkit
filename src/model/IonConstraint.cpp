#include "IonConstraint.h"
#include <string>


using namespace std;
/*************************
 * ION_CONSTRAINT methods
 *************************/

/**
 *Initializes an IonConstraint object, called from constructor 
 */
void IonConstraint::init(){

  use_neutral_losses_ = false;
  
  
  int modification_idx = 0;

  // initialize all modifications count to 0
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    modifications_[modification_idx] = 0;
  }

  mass_type_ = MONO;
  max_charge_ = 0;
  ion_type_ = (ION_TYPE_T)0;
  precursor_ion_ = false;
  min_charge_ = 0;
  exact_modifications_ = false;
  pointer_count_ = 0;

}

/**
 *\returns an empty heap allocated ion_constraint
 */
IonConstraint::IonConstraint() {

  init();
}

/**
 * modification, all modifications 0
 * add more modifications as needed using the set_ion_constraint_modification
 *\returns a new heap allocated ion_constraint
 */
IonConstraint::IonConstraint(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< max charge of the ions <= the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  bool precursor_ion  ///< should include precursor ion?
  )
{

  init();
  use_neutral_losses_ = false;
 
  // set all fields of constraint
  mass_type_ = mass_type;

  string charge_str = get_string_parameter_pointer("max-ion-charge");

  max_charge_ = max(1, max_charge - 1);

  if (charge_str != string("peptide")) {
    int charge_val;
    bool success = from_string(charge_val, charge_str);
    if (success) {
      max_charge_ = min(charge_val, max_charge_);
    } else {
      carp_once(CARP_WARNING, "Charge is not valid:%s", charge_str.c_str());
    }
  }

  min_charge_ = 0;
  exact_modifications_ = false;
  ion_type_ = ion_type;
  precursor_ion_ = precursor_ion;
  pointer_count_ = 1;

}

/**
 * \brief Create a new ion constraint based on the score type and the
 * charge of the peptide to be modeled.  Uses other
 * new_ion_constraint_ methods for some types.
 *
 * \returns A newly allocated ion constraint.
 */
IonConstraint* IonConstraint::newIonConstraintSmart(
  SCORER_TYPE_T score_type,
  int charge
){

  IonConstraint* new_constraint = NULL;

  switch(score_type){
  case SP:
    new_constraint = newIonConstraintSequestSp(charge);
    break;
  case XCORR:
    new_constraint = newIonConstraintSequestXcorr(charge);
    break;
  default:
    // use default type for others
    new_constraint = new IonConstraint(
      get_mass_type_parameter("fragment-mass"),
      charge,
      get_ion_type_parameter("primary-ions"),
       get_boolean_parameter("precursor-ions")); 
    break;
  }
  return new_constraint;
}

/**
 * modification, sets all fields for gmtk settings
 *\returns a new heap allocated ion_constraint
 */
IonConstraint* IonConstraint::newIonConstraintGmtk(
  int charge 
  )
{
  IonConstraint* constraint = NULL;

  int max_charge = 1;
  if(charge > 1){
    max_charge = charge;
  }  
  constraint = new IonConstraint(MONO, max_charge, ALL_ION, false);

  // set all modifications count for gmtk
  constraint->use_neutral_losses_ = true;
  constraint->min_charge_ = 1;
  constraint->modifications_[NH3] = -1;
  constraint->modifications_[H2O] = -1;
  constraint->modifications_[ISOTOPE] = 0;
  constraint->modifications_[FLANK] = 0;

  return constraint;
}


/**
 * modification, sets all fields for sequest settings
 *\returns a new heap allocated ion_constraint
 */
IonConstraint* IonConstraint::newIonConstraintSequest(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions. 
                  ///< cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  bool precursor_ion ///< should include precursor ion?
  )
{
  IonConstraint* constraint = new IonConstraint(mass_type, max_charge, ion_type, precursor_ion);

  // set                                                     
  constraint->use_neutral_losses_ = true;
  
  // set all modifications count for sequest
  constraint->modifications_[NH3] = 1;
  constraint->modifications_[H2O] = 1;
  constraint->modifications_[ISOTOPE] = 1;
  constraint->modifications_[FLANK] = 1;
  
  return constraint;
}

/**
 * modification, sets all fields for sequest Sp scoring settings
 * make B, Y type ions
 *\returns a new heap allocated ion_constraint
 */
IonConstraint* IonConstraint::newIonConstraintSequestSp(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  IonConstraint* constraint = NULL;
  constraint = new IonConstraint(MONO, charge, BY_ION, false);

  // set                                                     
  constraint->use_neutral_losses_ = true;
  
  // set all modifications count for sequest
  constraint->modifications_[NH3] = 0;
  constraint->modifications_[H2O] = 0;
  constraint->modifications_[ISOTOPE] = 0;
  constraint->modifications_[FLANK] = 0;
  
  return constraint;
}


/**
 * modification, sets all fields for Sequest Xcorr scoring settings
 * make B, Y, A type ions
 *\returns a new heap allocated ion_constraint
 */
IonConstraint* IonConstraint::newIonConstraintSequestXcorr(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  IonConstraint* constraint = NULL;
  constraint = new IonConstraint(MONO, charge, BYA_ION, false);

  // set                                                     
  constraint->use_neutral_losses_ = true;
  
  // set all modifications count for sequest
  constraint->modifications_[NH3] = 0;// -1;
  constraint->modifications_[H2O] = 0;// -1;
  constraint->modifications_[ISOTOPE] = 0;// not sure
  constraint->modifications_[FLANK] = 0;
  
  return constraint;
}

/**
 * Frees an allocated ion_constraint object.
 */
void IonConstraint::free(IonConstraint* ion_constraint)
{
  ion_constraint->pointer_count_--;
  if (ion_constraint->pointer_count_ == 0) {
    delete ion_constraint;
  }
}


IonConstraint::~IonConstraint()
{
}


/**
 * copies ion_constraint pointer
 */
IonConstraint* IonConstraint::copy(
    IonConstraint* ion_constraint
    ){
  
  ion_constraint->pointer_count_++;
  return ion_constraint;
}


/**
 * copies ion_constraint object from src to dest
 * must pass in a memory allocated ION_CONSTRAINT_T dest
 */
void IonConstraint::copy(
  IonConstraint* src,///< ion_constraint to copy from -in
  IonConstraint* dest///< ion_constraint to copy to -out
  )
{
  int modification_idx = 0;
  dest->use_neutral_losses_ = src->use_neutral_losses_;

  // if use natural loss, copy
  if(src->use_neutral_losses_){
    // iterate over all modifications a update new constraint
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      dest->modifications_[modification_idx] = src->modifications_[modification_idx];
    }
  }
  
  dest->mass_type_ = src->mass_type_;
  dest->max_charge_ = src->max_charge_;
  dest->ion_type_ = src->ion_type_;
  dest->precursor_ion_ = src->precursor_ion_;
}

/** FIXME!!!! double check
 * Determines if a ion satisfies a ion_constraint.
 * \returns true if the constraint is satisified. false if not.
 */
bool IonConstraint::isSatisfied(
   Ion* ion ///< query ion -in
   )
{
  int* counts = NULL;

  // TODO Fix
  bool return_val = true;
  // print_ion(ion, stderr);
  // fprintf(stderr, "%i->%i\n", ion_constraint->min_charge, ion_constraint->max_charge);
  // check ion type
  ION_TYPE_T ion_type = ion->getType();
  if(
     !(ion_type == ion_type_)
      
     &&
     
     !((ion_type_ == BY_ION) && (ion_type == B_ION || ion_type == Y_ION)) 
     
     &&

     !((ion_type_ == BYA_ION) 
          && 
       (ion_type == B_ION || ion_type == Y_ION || ion_type == A_ION)) 
     
     &&
     
     !(ion_type_ == ALL_ION)

     ){
     
    // precursor ion?
    if(!(precursor_ion_ && ion_type == P_ION)){
      return_val = false;
    }
  }
  
  // check charge
  if(ion->getCharge() > max_charge_){
    return_val = false;
  }
  
  if(ion->getCharge() < min_charge_){
    return_val = false;
  }

  // check modifications
  counts = ion->getModificationCounts();
  int mod_idx;
  for(mod_idx=0; mod_idx < MAX_MODIFICATIONS; ++mod_idx){
    if(modifications_[mod_idx] >= 0){
      if(counts[mod_idx] > modifications_[mod_idx]){
        return_val = false;
        break;
      }
    }
    else{
      if(counts[mod_idx] < modifications_[mod_idx]){
        return_val = false;
        break;
      }
    }
    if (exact_modifications_){
      if(counts[mod_idx] != modifications_[mod_idx]){
        return_val = false; 
        break;
      }
    }
  }
  
  // Add more checks here as more contraints are added

  // fprintf(stderr, "r = %i\n", return_val);
  return return_val;
}

/**
 * \returns ION_TYPE for this constraint
 */
ION_TYPE_T IonConstraint::getIonType() {
  return ion_type_;
}



/**
 * sets the modification count
 * can only add isotopes
 */
void IonConstraint::setModification(
  ION_MODIFICATION_T mod_type, ///< ion modification type -in
  int count  ///< the count of the modification -in  
  )
{
  // set modification count, can only add for isotope
  if(mod_type != ISOTOPE){
    modifications_[mod_type] = count;
  }
  else{
    modifications_[mod_type] = abs(count);
  }

  // set neutral loss to true if needed
  if(!use_neutral_losses_){
    use_neutral_losses_ = true;
  }
}

/**
 * sets the exact modification boolean 
 */
void IonConstraint::setExactness(
  bool exactness ///< whether to use exact mods or not -in
  ){
  exact_modifications_ = exactness;
  if (exactness){
    min_charge_ = max_charge_;
  }
}
 
/**
 * gets the modification count for specific mod_type
 */
int IonConstraint::getModification(
  ION_MODIFICATION_T mod_type ///< ion modification type -in
  )
{
  return modifications_[mod_type];
}

/**
 * \returns the modifications array
 */
int* IonConstraint::getModifications() {
  return modifications_;
}


/**
 * gets the mass type of the ion_constraint
 */
MASS_TYPE_T IonConstraint::getMassType()
{
  return mass_type_;
}

/**
 * get the maximum charge of the IonConstraint
 */
int IonConstraint::getMaxCharge() 
{
  return max_charge_;
}

/**
 *\returns precuror_ion_
 */
bool IonConstraint::getPrecursorIon() {
  return precursor_ion_;
}

/**
 * get the neutral loss field of the ion constraint.
 */
bool IonConstraint::getUseNeutralLosses()
{
  return use_neutral_losses_;
}

