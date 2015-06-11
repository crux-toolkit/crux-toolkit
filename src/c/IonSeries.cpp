/*************************************************************************//**
 * \file ion_series.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 21 Sep 2006
 * \brief code to support working with a series of ions
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "Ion.h"
#include "IonSeries.h"
#include "utils.h"
#include "crux-utils.h"
#include "parameter.h"
#include "Peptide.h"
#include "mass.h"

#include "IonFilteredIterator.h"
#include "Spectrum.h"

using namespace Crux;

static const int BINARY_GMTK = 1;
static const int PRINT_NULL_IONS = 1;
static const int MIN_FRAMES = 3;



/**
 * \struct loss_limit
 * \brief An object that specifies the max amount of neutral loss
 * possible at a given cleavage index. 
 * All numbers are for forward ions(A,B,C) subtract from total to get
 * reverse limit.
 */
struct loss_limit{
  int nh3; ///< the limit to how many NH3 may be lost
  int h2o; ///< the limit to how many H2O may be lost
  // add more if needed for other neutral loss
  // If change this struct, must also modify update_ion_series method
};


/**
 * Initializes an (empty) ion_series object, only called by constructors.
 */
void IonSeries::init() {

  peptide_ = NULL;
  modified_aa_seq_ = NULL;
  peptide_mass_ = 0;
  charge_ = 0;
  constraint_ = NULL;
  is_predicted_ = false;
  
  loss_limit_ = NULL;
  peptide_length_ = 0;
}

/**
 * \returns An (empty) ion_series object.
 */
IonSeries::IonSeries() {
  init();
}

/**
 * \brief Creates an ion series for a specific peptide, charge state,
 * and constraint without acutally predicting the ions.
 *
 * Creates a copy of the peptide sequence and calcualtes the peptide
 * mass.
 * Use this method to create ion_series only when few are needed
 * because the memory allocation process is expensive.  Alternatively,
 * use "new_ion_series_generic" and "update_ion_series" in combination
 * to reuse an ion_seires object.
 * \returns A newly allocated ion_series object from the given
 * peptide, charge, and constraint.
 */
IonSeries::IonSeries(
  const char* peptide, ///< The peptide for this ion series. -in
  int charge, ///< The charge for this ion series -in
  IonConstraint* constraint ///< constraints which these ions obey.
  )
{

  init();
  // copy the peptide sequence
  peptide_ = my_copy_string(peptide);
  convert_to_mod_aa_seq(peptide_, &(modified_aa_seq_));
  peptide_mass_ = Peptide::calcSequenceMass(peptide_, MONO);
  charge_ = charge;
  constraint_ = constraint;
  peptide_length_ = strlen(peptide_);
  
  // create the loss limit array
  loss_limit_ = 
    (LOSS_LIMIT_T*)mycalloc(peptide_length_, sizeof(LOSS_LIMIT_T));
}

/**
 * \brief Creates an ion_series object without a specific peptide.
 * Peptide details are added with the "update_ion_series" method.
 * \returns A newly allocated ion_series object that must be updated
 * for each peptide instance.
 */
IonSeries::IonSeries(
  IonConstraint* constraint, ///< constraints which these ions obey.
  int charge ///< The charge for this ion series -in
  ) {

  init();
  constraint_ = constraint;
  charge_ = charge;
  // use max peptide len so loss_limit array can be used for any peptide
  loss_limit_ = 
    (LOSS_LIMIT_T*)mycalloc(get_int_parameter("max-length"), 
                            sizeof(LOSS_LIMIT_T));
}

/**
 * \brief Updates an ion_series to a specific instance of a peptide
 * sequence. If the ion_series has already generated ions, they will
 * be free-ed. A copy of the peptide sequence is made and all other
 * variables (except ion_constraint) are also updated for the new
 * peptide sequence. 
 */
void IonSeries::update(
  char* peptide, ///< The peptide sequence with no mod characters. -in
  MODIFIED_AA_T* mod_seq ///< modified version of char* sequence -in
  ) 
{
  int ion_type_idx = 0;

  // Initialize the ion_series object for the new peptide sequence
  
  // free old peptide sequence
  if(peptide_){
    free(peptide_);
  }
  if(modified_aa_seq_){
    free(modified_aa_seq_);
  }
  if(loss_limit_){
    free(loss_limit_);
  }
  
  // iterate over all ions, and free them
  for (unsigned int idx=0;idx<ions_.size();idx++) {
    delete ions_[idx];
  }
  ions_.clear();
  
  // initialize all specific_ions back to 0
  int ion_idx;

  for(ion_idx=0; ion_type_idx < MAX_NUM_ION_TYPE; ++ion_type_idx){
    specific_ions_[ion_type_idx].clear();
  }

  is_predicted_ = false;
  
  // set ion_series for new instance of peptide
  
  // copy the peptide sequence
  peptide_ = my_copy_string(peptide);
  peptide_length_ = strlen(peptide);
  modified_aa_seq_ = copy_mod_aa_seq(mod_seq, peptide_length_);
  
  // Initialize the loss limit array for the new peptide
  loss_limit_ = 
    (LOSS_LIMIT_T*)mycalloc(peptide_length_, sizeof(LOSS_LIMIT_T));
  for(ion_idx =0; ion_idx < peptide_length_; ++ion_idx){
    loss_limit_->h2o = 0;
    loss_limit_->nh3 = 0;
    // add more initialize count if more added
  }
}


/**
 * Frees an allocated ion_series object.
 */
IonSeries::~IonSeries()
{
  if(peptide_){
    free(peptide_);
  }
  if(modified_aa_seq_){
    free(modified_aa_seq_);
  }
  if(loss_limit_){
    free(loss_limit_);
  }
  // free constraint?

  // iterate over all ions, and free them

  for (unsigned int idx=0;idx<ions_.size();idx++) {
    Ion::freeIon(ions_[idx]);
  }
  ions_.clear();

}


/**
 * Iterator access
 */
IonIterator IonSeries::begin() {
  return ions_.begin();
}

IonIterator IonSeries::end() {
  return ions_.end();
}

/**
 * Prints a ion_series object to file.
 */
void IonSeries::print(
  FILE* file ///< file for output -out
  )
{
  // check if the ions has been already predicted
  if(!is_predicted_){
    carp(CARP_ERROR, "ion series has not predicted ions");
    return;
  }
  
  // print header
  fprintf(file, "m/z\tmass\tcharge\tion-series\tpeptide-bond-index\tnh3\th2o\tisotope\tflank\n");
  
  
  // print each ion in the ion series
  unsigned int ion_idx;
  for(ion_idx=0; ion_idx < ions_.size(); ++ion_idx){
    ions_[ion_idx]->print(file);
  }
}

/**
 * Prints a ion_series object to file, in GMTK single-ion format.
 */
void IonSeries::printSingleGmtk(
  IonConstraint* ion_constraint, ///< ion_constraint to obey -in 
  FILE* file,                       ///< file output
  int sentence_idx){

  // create the filtered iterator that will select among the ions
  IonFilteredIterator ion_iterator(this, ion_constraint);
  
  // foreach ion in ion iterator, add matched observed peak intensity
  Ion* ion;
  int frame_idx = 0;
  while(ion_iterator.hasNext()) {
    ion = ion_iterator.next();
    
#ifdef BINARY_GMTK
    ion->printGmtkSingleBinary(file, sentence_idx, frame_idx);
#else
    ion->printGmtkSingle(file);
    sentence_idx++; // hack to avoid error for not using sentence_idx
#endif
    frame_idx++;
  }
  
  // print a null ion if there are none in this ion series
#ifdef PRINT_NULL_IONS
  for (; frame_idx < MIN_FRAMES; frame_idx++){
#ifdef BINARY_GMTK
    print_null_ion_gmtk_single_binary(file, sentence_idx, frame_idx);
#else
    print_null_ion_gmtk_single(file);
    sentence_idx++; // hack to avoid error for not using sentence_idx
#endif
  }
#endif

}

/**
 * Prints a ion_series object to file, in GMTK paired-ion format.
 */
void IonSeries::printPairedGmtk(
  IonConstraint* first_ion_constraint, ///< ion_constraint to obey -in 
  IonConstraint* second_ion_constraint, ///< ion_constraint to obey -in 
  FILE* file, ///< file output
  int sentence_idx
  )
{
  
  // create the filtered iterator that will select among the ions
  IonFilteredIterator ion_iterator(this, first_ion_constraint);
  
  // foreach ion in ion iterator, add matched observed peak intensity
  int frame_idx = 0;
  while(ion_iterator.hasNext()){
    Ion* first_ion = ion_iterator.next();
    int cleavage_idx = first_ion->getCleavageIdx();
    Ion* second_ion = this->getIon(second_ion_constraint, cleavage_idx);
    if ( (first_ion == NULL) || (second_ion == NULL) ){
      continue;
      }
    Ion::printGmtkPairedBinary(
                                 first_ion, 
                                 second_ion, 
                                 file,
                                 sentence_idx,
                                 frame_idx++);
  }
  
#ifdef PRINT_NULL_IONS
  for (; frame_idx < MIN_FRAMES; frame_idx++){
    print_null_ion_gmtk_paired_binary(file, sentence_idx, frame_idx);
  }
#endif

}


/**
 * \brief Find instances of amino acid which can incur neutral
 * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
 * Set the count of those observed so far for each cleavage index.
 * If no instance of amino acid, the count is assigned to 0
 * The information is used to determine if how many nh3 or h2o neutral
 * losses are possible. 
 */
void IonSeries::scanForAAForNeutralLoss()
{
  int peptide_length = peptide_length_;
  char* sequence = peptide_;

  // make sure loss_limit array is the right size
  if (peptide_length_ != strlen(sequence)){
    if (loss_limit_){
      free(loss_limit_);
    }
    loss_limit_ = 
      (LOSS_LIMIT_T*)mycalloc(peptide_length_, sizeof(LOSS_LIMIT_T));
  }

  int h2o_aa = 0;
  int nh3_aa = 0;

  // search for the first instance of the amino acids
  int cleavage_idx;
  for(cleavage_idx=0; cleavage_idx < peptide_length; ++cleavage_idx){
    LOSS_LIMIT_T* loss_limit_count = &loss_limit_[cleavage_idx];
    // is the AA  (S|T|E|D) ?
    if(sequence[cleavage_idx] == 'S' ||
       sequence[cleavage_idx] == 'T' ||
       sequence[cleavage_idx] == 'E' ||
       sequence[cleavage_idx] == 'D' )
      {
        loss_limit_count->h2o = ++h2o_aa;
        loss_limit_count->nh3 = nh3_aa;
      }
    // is the AA  (R|K|Q|N) ?
    else if(sequence[cleavage_idx] == 'R' ||
            sequence[cleavage_idx] == 'K' ||
            sequence[cleavage_idx] == 'Q' ||
            sequence[cleavage_idx] == 'N' )
      {
        loss_limit_count->nh3 = ++nh3_aa;
        loss_limit_count->h2o = h2o_aa;
      }
    else{
      loss_limit_count->h2o = h2o_aa;
      loss_limit_count->nh3 = nh3_aa;
    }
  }

}

/**
 * \brief Creates an array in which element i is the sum of the masses
 * of amino acids 0 to (i-1).  At i=0 is stored the length of the
 * peptide.  
 * \returns an array of ion masses for all sub sequences
 */
FLOAT_T* IonSeries::createIonMassMatrix(
  //char* peptide, ///< The peptide for this ion series. -in
  MODIFIED_AA_T* modified_seq, ///< the sequence
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int peptide_length ///< the length of the peptide
  )
{
  if( modified_seq == NULL ){
  //if( peptide == NULL ){
    carp(CARP_ERROR, "Cannot create mass matrix from NULL seqence");
    return NULL;
  }

  FLOAT_T* mass_matrix = (FLOAT_T*)mymalloc(sizeof(FLOAT_T)*(peptide_length+1));
  
  // at index 0, the length of the peptide is stored
  mass_matrix[0] = peptide_length;

  // add up AA masses
  int ion_idx = 1;
  // initialize first to be mass of c-term amino acid
  // mass_matrix[ion_idx] = get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  mass_matrix[ion_idx] = get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
  //++ion_idx;
  //for(; ion_idx <= peptide_length; ++ion_idx){
  for(ion_idx = 2; ion_idx <= peptide_length; ++ion_idx){
    mass_matrix[ion_idx] = mass_matrix[ion_idx-1] + 
      get_mass_mod_amino_acid(modified_seq[ion_idx-1], mass_type);
      //get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  }
  // DEBUGGING
  /*
  fprintf(stderr, "seq:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%s", modified_aa_to_string(modified_seq[ion_idx]));
  }
  fprintf(stderr, "\nmas:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", get_mass_mod_amino_acid(modified_seq[ion_idx], MONO));
  }
  fprintf(stderr, "\nsum:");
  for(ion_idx = 0; ion_idx < peptide_length; ++ion_idx){
    fprintf(stderr, "\t%.2f", mass_matrix[ion_idx+1]);
  }
  fprintf(stderr, "\n");
  */
  return mass_matrix;
}

/**
 * user must ensure that there is enough space for this ion
 * adds ion to ion_series' master ion_array and if B|Y ion to the specific ion_array
 */
void IonSeries::addIon(
  Ion* ion ///< ion to add -in
  )
{
  // add ion to ion series
  ions_.push_back(ion);
  
  // increment the pointer
  ion->incrementPointerCount();

  // add a pointer of ion to the specific ion_type array
  specific_ions_[ion->getType()].push_back(ion);
}

/**
 * helper function: add_ions
 * add all the ions to ion_series up to the max charge
 *\returns true if successfully adds all ions, else false
 */
bool IonSeries::addIonsByCharge(
  FLOAT_T mass, ///< the base mass of the ion to add
  int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
  ION_TYPE_T ion_type ///< the ion type of the ions to be added
  )
{
  IonConstraint* constraint = constraint_;
  int charge_idx = 1;
  Ion* ion = NULL;
  int max_charge;
  
  // set the max charge, the maximum cannot exceed the precursor ion's charge
  if(constraint->getMaxCharge() > charge_){
    max_charge = charge_;
  }
  else{
    max_charge = constraint->getMaxCharge();
  }

  // iterate over all different charge
  for(; charge_idx <= max_charge; ++charge_idx){
    // create ion
    ion = new Ion(ion_type, 
                  cleavage_idx, 
                  charge_idx, 
                  peptide_, 
                  constraint->getMassType(), 
                  mass); 
    // add ion to ion series
    this->addIon(ion);
  }
  
  return true;
}

/**
 * Creates all the ions with no modifications up to the max charge
 * Adds each ion to ion_series
 *\returns true if successfully generates all the ions, else false
 */
bool IonSeries::generateIonsNoModification(
  FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
  )
{
  if( mass_matrix == NULL ){
    carp(CARP_ERROR,
         "Cannot generate ions from NULL mass matrix");
    return false;
  }
  int cleavage_idx = 1;
  IonConstraint* constraint = constraint_;
  FLOAT_T mass = 0;

  // get peptide length
  int peptide_length = (int)mass_matrix[0];

  // iterate over all cleavage index
  for(; cleavage_idx < peptide_length; ++cleavage_idx){
    
    // add A ion
    if(constraint->getIonType() == A_ION 
       || constraint->getIonType() == BYA_ION 
       || constraint->getIonType() == ALL_ION){

      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->getMassType() == MONO){
        mass -= MASS_CO_MONO; 
      }
      else{ // average
        mass -= MASS_CO_AVERAGE; 
      }

      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, A_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for A ion");
      return false;
      }
    }
    
    // add B ion
    if(constraint->getIonType() == ALL_ION 
       || constraint->getIonType() == BY_ION
       || constraint->getIonType() == BYA_ION
       || constraint->getIonType() == B_ION){
      
      // set mass
      mass = mass_matrix[cleavage_idx];

      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, B_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for B ion");
        return false;
      }

    }
    
    // add C ion
    if(constraint->getIonType() == C_ION || constraint->getIonType() == ALL_ION){
      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->getMassType() == MONO){
        mass += MASS_NH3_MONO; 
      }
      else{ // average
        mass += MASS_NH3_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, C_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for C ion");
        return false;
      }
    }
    
    // add X ion
    if(constraint->getIonType() == X_ION || constraint->getIonType() == ALL_ION){
      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];

      if(constraint->getMassType() == MONO){
        mass += MASS_CO_MONO + MASS_H2O_MONO;     
      }
      else{ // average
        mass += MASS_CO_AVERAGE + MASS_H2O_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, X_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for X ion");
        return false;
      }
    }
    
    // add Y ion
    if(constraint->getIonType() == ALL_ION || 
       constraint->getIonType() == BY_ION ||
       constraint->getIonType() == BYA_ION ||
       constraint->getIonType() == Y_ION){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint->getMassType() == MONO){
        mass += MASS_H2O_MONO; 
      }
      else{ // average
        mass += MASS_H2O_AVERAGE;
      }

      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, Y_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Y ion");
        return false;
      }
      
    }
    
    // add Z ion
    if(constraint->getIonType() == Z_ION ||
       constraint->getIonType() == ALL_ION ){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint->getMassType() == MONO){
        mass = mass - MASS_NH3_MONO + MASS_H2O_MONO; 
      }
      else{ // average
        mass = mass - MASS_NH3_AVERAGE + MASS_H2O_AVERAGE;
      }
      
      // add ions up to max charge
      if(!this->addIonsByCharge(mass, cleavage_idx, Z_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Z ion");
        return false;
      }
    }
  }

  // add P ion(precursor ion)?
  if(constraint_->getPrecursorIon()){
    
    // set mass 
    mass = mass_matrix[(int)mass_matrix[0]];
    
    // mass type
    if(constraint->getMassType() == MONO){
      mass += MASS_H2O_MONO; 
    }
    else{ // average
      mass += MASS_H2O_AVERAGE;
    }
    
    // add ions up to max charge
    if(!addIonsByCharge(mass, (int)mass_matrix[0], P_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge P ion");
      return false;
    }
  }

  return true;
}

/**
 * The modification depends on the loss/add && if the ion contains RKQN or STED
 * The number of losses possible cannot exceed the number of RKQN or STED in the ion
 * The loss_limit array in the ion_series must be populated prior to this method call
 *\returns true if the ion can lose the mod_type modification, else false
 */
bool IonSeries::canIonGenerateModification(
  Ion* ion, ///< the ion to check if can lose nh3 -in
  ION_MODIFICATION_T mod_type, ///< generate ions of this modification_type -in/out
  int increment  ///< the add/loss of the modification
  )
{

  //Make sure that an ion can have 0 or 1 modifications.
  if (ion->getTotalModificationCount() == 1) {
    return false;
  }
  

  int cleavage_idx = ion->getCleavageIdx();

  // check for NH3 modification
  if(mod_type == NH3){
    // adding is ok
    if(increment >= 0){
      return true;
    }
    
    // is forward ion_type(ABC)?
    if(ion->isForwardType()){
      // does this ion contain enough counts of the RKQN
      if(-increment >  (&loss_limit_[cleavage_idx-1])->nh3){
        return false;
      }
      return true;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the RKQN
      if(cleavage_idx == peptide_length_){
        if(-increment > (&loss_limit_[peptide_length_-1])->nh3){
          return false;
        }
      }
      else if(-increment >  
              ((&loss_limit_[peptide_length_-1])->nh3 - 
               (&loss_limit_[peptide_length_ - cleavage_idx - 1])->nh3)){
        return false;
      }
      return true;
    }
  }

  // check for H2O modification
  if(mod_type == H2O){
    // adding is ok
    if(increment >= 0){
      return true;
    }
    
    // is forward ion_type(ABC)?
    if(ion->isForwardType()){
      // does this ion contain enough counts of the STED
      if(-increment >  (&loss_limit_[cleavage_idx-1])->h2o){
        return false;
      }
      return true;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the STED
      if(cleavage_idx == peptide_length_){
        if(-increment > (&loss_limit_[peptide_length_-1])->h2o){
          return false;
        }
      }
      else if(-increment >  
              ((&loss_limit_[peptide_length_-1])->h2o - 
               (&loss_limit_[peptide_length_ - cleavage_idx - 1])->h2o)){
        return false;
      }
      return true;
    }
  }
  
  // check for ISOTOPE modification
  else if(mod_type == ISOTOPE){
    // adding is ok
    if(increment >= 0 && ion->getType() != P_ION){
      return true;
    }
    // add more constraint if needed
  }
  // check for FLANK modification
  else if(mod_type == FLANK){
    // only add flanking ions to type B,Y ions
    if(ion->getType() == B_ION || ion->getType() == Y_ION){
      // only add ions with no modifications
      if(!ion->isModified()){
        return true;
      }
    }
  }
  
  return false;
}

/**
 * creates all the ions with specific modifications up to the max charge
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns true if successfully generates all the ions with modifications, else false
 */
bool IonSeries::generateIons(
  ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
  )
{
  int ion_idx = 0;
  int total_ion = ions_.size();
  Ion* working_ion = NULL;
  Ion* new_ion = NULL;
  int* modifications = constraint_->getModifications();

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // reset modification increment, if mod_type loss
  if(modifications[mod_type] < 0){
    type_increment = -1;
  }
  
  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ions_[ion_idx];
    
    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!(canIonGenerateModification(working_ion, mod_type, type_increment))){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; abs(type_idx) <= abs(modifications[mod_type]); ){
      // copy the src ion, into new ion
      new_ion = new Ion();
      Ion::copy(working_ion, new_ion, working_ion->getPeptideSequence());
      
      // add the modification to the new ion
      new_ion->addModification(mod_type, type_idx, constraint_->getMassType());
      
      // add ion to ion_series
      addIon(new_ion);
     
      // can this ion generate a mod_type modification for the next count of modification?, 
      if(!(this->canIonGenerateModification(working_ion, mod_type, 
              (type_idx += type_increment)))){
        break;
      }
    }
  }
  return true;
}

/**
 * creates all the flanking ions up to the max charge
 * can only create flanking ions that are B|Y ions and don't have modification
 * assumes the ions with no modification all are at the begining of the ion[] in ion_series
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns true if successfully generates all the ions with modifications, else false
 */
bool IonSeries::generateIonsFlank()
{
  int ion_idx = 0;
  int total_ion = ions_.size();
  Ion* working_ion = NULL;
  Ion* new_ion = NULL;
  Ion* new_ion2 = NULL;
  int* modifications = constraint_->getModifications();

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ions_[ion_idx];
    
    // no more ions that are not modified, thus done
    if(working_ion->getSingleModificationCount(NH3) != 0){
      break;
    }

    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!canIonGenerateModification(working_ion, FLANK, type_increment)){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; type_idx <= modifications[FLANK]; type_idx += type_increment){
      // copy the src ion, into new ion
      new_ion = new Ion();
      new_ion2 = new Ion();
      Ion::copy(working_ion, new_ion, working_ion->getPeptideSequence());
      Ion::copy(working_ion, new_ion2, working_ion->getPeptideSequence());
      
      // add the modification to the new ion
      new_ion->addModification(FLANK, type_idx, constraint_->getMassType());
      new_ion2->addModification(FLANK, -type_idx, constraint_->getMassType());

      // add ion to ion_series
      addIon(new_ion);
      addIon(new_ion2);
    }
  }
  return true;
}

/**
 * \brief The engine of ion series. Predicts all the ions from the
 * peptide that meet the ion constraint. All predicted ions are stored
 * in the ion_series as ion objects. 
 */
void IonSeries::predictIons()
{
  if(is_predicted_){
    carp(CARP_WARNING, "The ion series has already been predicted and added");
    return;
  }

  IonConstraint* constraint = constraint_;
  
  // create a mass matrix
  FLOAT_T* mass_matrix = 
    createIonMassMatrix(modified_aa_seq_, constraint->getMassType(), peptide_length_);  
  /*
  printf("cumulative mass sum is:\n");
  int idx = 0;
  for(idx = 0; idx < mass_matrix[0]; idx++){
    printf("%i\t%f\n", idx, mass_matrix[idx]);
  }
  */
  // scan for the first and last  (S, T, E, D) and (R, K, Q, N), 
  // initialize to determine modification is ok.
  // the first, last of STED, RKQN are stored in ion_series.

  scanForAAForNeutralLoss();

  // generate ions without any modifications
  if(!generateIonsNoModification(mass_matrix)){
    carp(CARP_FATAL, "failed to generate ions, no modifications");
  }

  // create modification ions?
  if(constraint_->getUseNeutralLosses()){
    
    // generate ions with nh3 modification
    if(abs(constraint->getModifications()[NH3]) > 0){
      if(!generateIons(NH3)){
        carp(CARP_FATAL, "failed to generate ions, NH3 modifications");
      }
    }
    
    // generate ions with h2o modification
    if(abs(constraint->getModifications()[H2O]) > 0){
      if(!generateIons(H2O)){
        carp(CARP_FATAL, "failed to generate ions, H2O modifications");
      }
    }

    // generate ions with isotope modification
    if(constraint->getModifications()[ISOTOPE] > 0){
      if(!generateIons(ISOTOPE)){
        carp(CARP_FATAL, "failed to generate ions, ISOTOPE modifications");
      }
    }

    // generate ions with flank modification
    if(constraint->getModifications()[FLANK] > 0){
      if(!generateIonsFlank()){
        carp(CARP_FATAL, "failed to generate ions, FLANK modifications");
      }
    }
    
    // add more modifications here

  }
  
  // ion series now been predicted
  is_predicted_ = true;
  // free mass matrix
  free(mass_matrix);

}

/**
 * Copies ion_series object from src to dest.
 *  must pass in a memory allocated ION_SERIES_T* dest
 * does not copy the loss_limit.
 */
void IonSeries::copy(
  IonSeries* src,///< ion to copy from -in
  IonSeries* dest///< ion to copy to -out
  )
{
  Ion* src_ion = NULL;
  Ion* dest_ion = NULL;
  
  dest->peptide_ = my_copy_string(src->peptide_);
  dest->charge_ = src->charge_;
  dest->peptide_length_ = src->peptide_length_;
  //mod seq???

  // add copy of pointer ion constraint
  dest->constraint_ = src->constraint_;

  // add copy ion, add ion_filtered_iterator

  for (IonIterator iterator = src->begin();
    iterator != src->end();
    ++iterator) {

    src_ion = *iterator;
    // add ion
    dest_ion = new Ion();
    Ion::copy(src_ion, dest_ion, dest->peptide_);
    dest->addIon(dest_ion);
  }

  dest->is_predicted_ = true;
}

/**
 * remove an ion from IonSeries, does not free ion.
 */
void IonSeries::removeIon(
  Ion* ion ///<ion to remove
  ) {


  IonIterator ion_iter;

  for (ion_iter = begin();
    ion_iter != end();
    ++ion_iter) {

    if (*ion_iter == ion) {
      ions_.erase(ion_iter);
      return;
    }

    

  }

  carp(CARP_ERROR, "Cannot find ion to delete:%i", ion);

}


/*************************************
 * ION_SERIES_T: get and set methods
 ************************************/

/**
 * \returns the ion that meets the criteria or NULL
 * TODO possibly reimplement if linear scan is too slow
 */
Ion* IonSeries::getIon(
  IonConstraint* ion_constraint,
  int cleavage_idx
  ){

  IonFilteredIterator ion_iterator(this, ion_constraint);
  Ion* ion = NULL;

  while(ion_iterator.hasNext()){
    ion = ion_iterator.next();
    if(ion->getCleavageIdx() == cleavage_idx){
      return ion;
    }
  }
  return NULL;
}

/**
 *\returns the peptide length of which the ions are made
 */
int IonSeries::getPeptideLength()
{
  return peptide_length_;
}

/**
 * User should not free the peptide sequence seperate from the ion_series
 *\returns a pointer to the original parent peptide sequence of the ion_series object
 */
char* IonSeries::getPeptide()
{
  return peptide_;
}

/**
 * copies in the peptide sequence to heap allocated sequence.
 * set the parent peptide sequence of the ion_series object
 */
void IonSeries::setPeptide(
  char* peptide///< the peptide sequence to set -in
  )
{
  // free previous sequence
  if(peptide_ != NULL){
    free(peptide_);
  }
  peptide_ = my_copy_string(peptide);
}

/**
 *\returns the charge of the ion_series object
 */
int IonSeries::getCharge()
{
  return charge_;
}

/**
 * set the charge of the ion_series object
 */
void IonSeries::setCharge(
  int charge///< the charge of the ion -in
  )
{
  charge_ = charge;
}

/**
 * get the is_predicted field of the ion_series object
 */
bool IonSeries::getIsPredicted() {
  return is_predicted_;
}

/**
 * set the is_predicted field of the ion_series object
 */
void IonSeries::setIsPredicted(
  bool is_predicted///< the is_predicted field -in
  ) {
  is_predicted_ = is_predicted;
}

/**
 * get the modified_aa_seq of the ion_series object
 */
MODIFIED_AA_T* IonSeries::getModifiedAASeq()
{
  return modified_aa_seq_;
}

/**
 *\returns the constraint of the ion_series object
 */
IonConstraint* IonSeries::getIonConstraint()
{
  return constraint_;
}

/**
 * set the of the ion_series object
 */
void IonSeries::setIonConstraint(
  IonConstraint* constraint///<  -in
  )
{
  constraint_ = constraint;
}

/**
 *\returns the total number of ions in the ion_series object
 */
int IonSeries::getNumIons()
{
  return ions_.size();
}

/**
 *\returns the total number of ion_type in the ion_series object
 */
int IonSeries::getNumIonsOneType(
  ION_TYPE_T ion_type ///< the type of ions -in
  )
{
  return specific_ions_[ion_type].size();
}

vector<Ion*>& IonSeries::getSpecificIons(
  ION_TYPE_T ion_type ///< the type of ions -in
  )
{
  return specific_ions_[ion_type];
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
