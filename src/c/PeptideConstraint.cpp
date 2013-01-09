/*************************************************************************//**
 * \file PeptideConstraint.cpp
 * \brief Object for holding the peptide constraint information.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "crux-utils.h"
#include "objects.h"
#include "mass.h"
#include "Peptide.h"
#include "Protein.h"
#include "carp.h"
#include "PeptideConstraint.h"

using namespace Crux;

/**
 * Allocates a new (empty) peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
void PeptideConstraint::init() {
  carp(CARP_DETAILED_DEBUG, "Initializing peptide constraint");
  enzyme_ = (ENZYME_T)0;
  digestion_ = (DIGEST_T)0;
  min_mass_ = 0;
  min_length_ = 0;
  max_length_ = 0;
  num_mis_cleavage_ = 0;
  mass_type_ = (MASS_TYPE_T)0;
  num_pointers_ = 1;

}

PeptideConstraint::PeptideConstraint() {
  init();
}

/**
 * Instantiates a new peptide_constraint object.
 * \returns An allocated PEPTIDE_CONSTRAINT_T object.
 */
PeptideConstraint::PeptideConstraint(

  ENZYME_T enzyme, 
  DIGEST_T digest,
  FLOAT_T min_mass, ///< the minimum mass -in
  FLOAT_T max_mass, ///< the maximum mass -in
  int min_length, ///< the minimum length of peptide -in
  int max_length,  ///< the maximum lenth of peptide(max limit = 255) -in
  int num_mis_cleavage, ///< The maximum mis cleavage of the peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  )
{
  // max length must be less or equal than 255 becuase of the unsigned char limit of 255
  if(max_length > 255){
    carp(CARP_FATAL, "ERROR: cannot set max length higer than 255");
  }
  
  init();

  setEnzyme(enzyme);
  setDigest(digest);
  setMinMass(min_mass);
  setMaxMass(max_mass);
  setMinLength(min_length);
  setMaxLength(max_length);
  setNumMisCleavage(num_mis_cleavage);
  setMassType(mass_type);
}

/**
 * \brief Create a new peptide constraint and populate its values
 * based on those in parameter.c 
 * \returns A newly allocated peptide constraint.
 */
PeptideConstraint* PeptideConstraint::newFromParameters() {

  PeptideConstraint* new_constraint = new PeptideConstraint();

  new_constraint->setEnzyme(get_enzyme_type_parameter("enzyme"));
  new_constraint->setDigest(get_digest_type_parameter("digestion"));
  new_constraint->setMinMass(get_double_parameter("min-mass"));
  new_constraint->setMaxMass(get_double_parameter("max-mass"));
  new_constraint->setMinLength(get_int_parameter("min-length"));
  new_constraint->setMaxLength(get_int_parameter("max-length"));
  // TODO : change this after missed cleavage is an integer parameter
  // rather than boolean.
  new_constraint->setNumMisCleavage(get_int_parameter("missed-cleavages"));
  new_constraint->setMassType(get_mass_type_parameter("isotopic-mass"));

  return new_constraint;
}


/**
 * Copy peptide pointer and increment pointer count
 */
PeptideConstraint* PeptideConstraint::copyPtr(
  PeptideConstraint* constraint
  ) {

  constraint->num_pointers_++;
  return constraint;
}

// FIXME check the association..as long as there is one tryptic parent then true
// num_miss_cleavage is not implemented..add if needed
/** 
 * Determines if a peptide satisfies a peptide_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
bool PeptideConstraint::isSatisfied(
  Peptide* peptide ///< the query peptide -in
  ) {

  return (peptide->getLength() <= getMaxLength() &&
     peptide->getLength() >= getMinLength() &&
     peptide->getPeptideMass() <= getMaxMass() &&
     peptide->getPeptideMass() >= getMinMass()
     );
}

/**
 * Frees an allocated peptide_constraint object.
 */
void PeptideConstraint::free(
  PeptideConstraint* peptide_constraint ///< object to free -in 
  ) {

  peptide_constraint->num_pointers_--;
  if (peptide_constraint->num_pointers_ <= 0) {
    carp(CARP_DETAILED_DEBUG, "Final free of peptide constraint");
    delete peptide_constraint;
  }

}

PeptideConstraint::~PeptideConstraint() {

}

/**
 * Sets the enzyme used for the in silicos digestion
 * of the protein sequence into peptides.
 */
void PeptideConstraint::setEnzyme(
  ENZYME_T enzyme
){

  enzyme_ = enzyme;
}

/**
 * \returns The enzyme for this peptide constraint.
 */
ENZYME_T PeptideConstraint::getEnzyme()
{
  return enzyme_;
}

/**
 * Sets the level of digestion for the peptide constraint.
 */
void PeptideConstraint::setDigest(
  DIGEST_T digest
  ){

  digestion_ = digest;
}

/**
 * \returns The level of digestion for the peptide constraint.
 */
DIGEST_T PeptideConstraint::getDigest() {

  return digestion_;
}


/**
 * sets the min mass of the peptide_constraint
 */
void PeptideConstraint::setMinMass(
  FLOAT_T min_mass  ///< the min mass of the peptide constraint - in
  )
{
  min_mass_ = min_mass;
}

/**
 * \returns the min mass of the peptide_constraint
 */
FLOAT_T PeptideConstraint::getMinMass() {

  return min_mass_;
}

/**
 * sets the max mass of the peptide_constraint
 */
void PeptideConstraint::setMaxMass(
  FLOAT_T max_mass  ///< the max mass of the peptide constraint - in
  ) {

  max_mass_ = max_mass;
}

/**
 * \returns the max mass of the peptide_constraint
 */
FLOAT_T PeptideConstraint::getMaxMass() {

  return max_mass_;
}

/**
 * sets the min length of the peptide_constraint
 */
void PeptideConstraint::setMinLength(
  int min_length  ///< the min length of the peptide constraint - in
  ) {

  min_length_ = min_length;
}

/**
 * \returns the min length of the peptide_constraint
 */
int PeptideConstraint::getMinLength() {

  return min_length_;
}

/**
 * sets the max length of the peptide_constraint
 * maximum limit 255
 */
void PeptideConstraint::setMaxLength(
  int max_length  ///< the max length of the peptide constraint - in
  ) {

  // check if maximum length is with in range <= 255
  if(max_length > 255){
    carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
  }
  
  max_length_ = max_length;
}

/**
 * \returns the max length of the peptide_constraint
 */
int PeptideConstraint::getMaxLength(
  ) {

  return max_length_;
}


/**
 * sets the num_mis_cleavage of the peptide_constraint
 */
void PeptideConstraint::setNumMisCleavage(
  int num_mis_cleavage ///< The maximum mis cleavage of the peptide -in
  ) {

  num_mis_cleavage_ = num_mis_cleavage;
}

/**
 * \returns the num_mis_cleavage of the peptide_constraint
 */
int PeptideConstraint::getNumMisCleavage() {

  return num_mis_cleavage_;
}

/**
 * sets the mass type of the peptide_constraint
 */
void PeptideConstraint::setMassType(
  MASS_TYPE_T mass_type ///< the mass_type for the constraint -in
  ) {

  mass_type_ = mass_type;
}

/**
 * \returns the mass type of the mass_constraint
 */
MASS_TYPE_T PeptideConstraint::getMassType() {

  return mass_type_;
}
