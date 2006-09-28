/*****************************************************************************
 * \file ion_series.c
 * AUTHOR: Chris Park
 * CREATE DATE: 21 Sep 2006
 * DESCRIPTION: code to support working with a series of ions
 * REVISION: $Revision: 1.6 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "ion.h"
#include "ion_series.h"
#include "utils.h"
#include "crux-utils.h"
#include "mass.h"

#define MAX_IONS 10000

/**
 * \struct ion_series
 * \brief An object to represent a series of ions, and organize them!
 * For which additional data structures will be created as needed 
 */
struct ion_series {
  char* peptide; ///< The peptide for this ion series
  int charge; ///< The charge state of the peptide for this ion series
  ION_CONSTRAINT_T* constraint; ///< The constraints which the ions in this series obey
  ION_T* ions[MAX_IONS]; ///< The ions in this series
  int num_ions; ///< the number of ions in this series
  BOOLEAN_T is_predicted; ///< has this ion_series been predicted ions already?
};

/**
 * \struct ion_constraint
 * \brief An object to represent the contraints which the ions in this
 * series obey.
 * CHRIS you can add BOOLEANS and other things as needed for implementation
 * of the predict-peptide-ions executable
 */
struct ion_constraint {
  BOOLEAN_T use_neutral_losses; ///< A boolean to determine if the ions series should include neutral losses
  int modifications[MAX_MODIFICATIONS]; ///< an array of to indicate which modifications to perform
  MASS_TYPE_T mass_type; ///< the mass_type to use MONO|AVERAGE
  int max_charge; ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type; ///< the ion types the peptide series should include
};

/**
 *\struct ion_iterator
 *\brief An object to iterate over all ion objects in the ion_series
 */
struct ion_iterator {
  ION_SERIES_T* ion_series; ///< the ion series that the ion we are iterating
  int ion_idx; ///< the current ion that is being returned 
};

/**
 *\struct ion_filtered_iterator
 *\brief An object to iterate over ion objects that meet constraint in the ion_series
 */
struct ion_filtered_iterator {
  ION_SERIES_T* ion_series; ///< the ion series that the ion we are iterating
  ION_CONSTRAINT_T* constraint; ///< The constraints which the ions in this series obey
  BOOLEAN_T has_next; ///< the boolean which the iterator has a next ion
  int ion_idx; ///< the current ion that is being returned 
  ION_T* ion; ///< the next ion to return when called upon
};

/**
 * \Returns An (empty) ion_series object.
 */
ION_SERIES_T* allocate_ion_series(void){
  ION_SERIES_T* ion_series = (ION_SERIES_T*)mycalloc(1, sizeof(ION_SERIES_T));
  ion_series->is_predicted = FALSE;
  
  return ion_series;
}

/**
 * copies in the peptide sequence
 * Instantiates a new ion_series object from a filename. 
 */
ION_SERIES_T* new_ion_series(
  char* peptide, ///< The peptide for this ion series. -in
  int charge, ///< The charge for this ion series -in
  ION_CONSTRAINT_T* constraint ///< The constraints which the ions in this series obey.
  )
{
  ION_SERIES_T* ion_series = allocate_ion_series();
  //copy the peptide sequence
  ion_series->peptide = my_copy_string(peptide);
  ion_series->charge = charge;
  ion_series->constraint = constraint;
  return ion_series;  
}

/**
 * Frees an allocated ion_series object.
 */
void free_ion_series(
  ION_SERIES_T* ion_series ///< the ion collection to free - in
  )
{
  free(ion_series->peptide);
  
  //iterate over all ions, and free them
  while(ion_series->num_ions > 0){
    free_ion(ion_series->ions[ion_series->num_ions-1]);
    --ion_series->num_ions;
  }
  
  free(ion_series);
}

/**
 * Prints a ion_series object to file.
 */
void print_ion_series(
  ION_SERIES_T* ion_series, ///< ion_series to print -in 
  FILE* file ///< file for output -out
  )
{
  //check if the ions has been already predicted
  if(!ion_series->is_predicted){
    carp(CARP_ERROR, "ion series has not predicted ions");
    return;
  }
  
  //print header
  fprintf(file, "m/z\tmass\tcharge\tion-series\tpeptide-bond-index\tNH3\tH2O\tISOTOPE\tFLANK\n");
  
  int ion_idx = 0;
  
  //print each ion in the ion series
  for(; ion_idx < ion_series->num_ions; ++ion_idx){
    print_ion(ion_series->ions[ion_idx], file);
  }
}

/**
 * index starts at 1, the length of the peptide is stored at index 0.
 * heap allocated mass matrix must be freed by user
 *\returns a ion mass matrix of all possible length
 */
float* create_ion_mass_matrix(
  char* peptide, ///< The peptide for this ion series. -in
  MASS_TYPE_T mass_type ///< the mass_type to use MONO|AVERAGE
  )
{
  int peptide_length = strlen(peptide);
  float* mass_matrix = (float*)mymalloc(sizeof(float)*(peptide_length+1));
  
  //at index 0, the length of the peptide is stored
  mass_matrix[0] = peptide_length;

  //add up AA masses
  int ion_idx = 1;
  mass_matrix[ion_idx] = get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  ++ion_idx;
  for(; ion_idx <= peptide_length; ++ion_idx){
    mass_matrix[ion_idx] = mass_matrix[ion_idx-1] + get_mass_amino_acid(peptide[ion_idx-1], mass_type);
  }
  return mass_matrix;
}

/**
 * helper function: add_ions
 * add all the ions to ion_series up to the max charge
 *\returns TRUE if successfully adds all ions, else FALSE
 */
BOOLEAN_T add_ions_by_charge(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in
  float mass, ///< the base mass of the ion to add
  int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
  int* modification_counts, ///< the different modifications required
  ION_TYPE_T ion_type ///< the ion type of the ions to be added
  )
{
  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  int charge_idx = 1;
  ION_T* ion = NULL;
  int max_charge;

  //check if there's enough space to add the new ions
  if(ion_series->num_ions + constraint->max_charge > MAX_IONS){
    carp(CARP_ERROR, "exceeds ion array size in ion_series");
    return FALSE;
  }
  
  //set the max charge, the maximum cannot exceed the precursor ion's charge
  if(constraint->max_charge > ion_series->charge){
    max_charge = ion_series->charge;
  }
  else{
    max_charge = constraint->max_charge;
  }

  //iterate over all different charge ////FIXME..make sure that max charge isn't too high!!!!!!!!!!
  for(; charge_idx <= max_charge; ++charge_idx){
    //create ion
    ion = new_modified_ion_with_mass(ion_type, cleavage_idx, charge_idx, ion_series->peptide, 
                                     constraint->mass_type, mass, modification_counts); 
    //add ion to ion series
    ion_series->ions[ion_series->num_ions++] = ion;

    //debug//
    /*
    ++TOTAL_IONS_ADDED;
    print_ion(ion, stdout);
    printf("total ions: %d\n", TOTAL_IONS_ADDED);
    */
  }

  return TRUE;
}

/**
 * helper function: predict_ions
 * adds ions to the ion series according to their type
 *\returns TRUE is successfully adds all the types of ions, ELSE flase
 */
BOOLEAN_T add_ions(         
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in/out
  int cleavage_idx, ///< the cleavage indx starting from the left -in
  float* mass_matrix, ///< the mass matrix that hold all ions base mass -in
  int* modification_counts ///< the different modifications required -in
  )
{
  float mass = 0;
  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  ION_TYPE_T ion_type = constraint->ion_type;

  //add A ion
  if(ion_type == A_ION){
    //set mass
    mass = mass_matrix[cleavage_idx];

    if(constraint->mass_type == MONO){
      mass -= MASS_CO_MONO; 
    }
    else{ //average
      mass -= MASS_CO_AVERAGE; 
    }

    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, A_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge for A ion");
      return FALSE;
    }
  }
  
  //add B ion
  if(ion_type == ALL_ION || ion_type == B_ION){
    //set mass
    mass = mass_matrix[cleavage_idx];

    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, B_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge for B ion");
      return FALSE;
    }
  }

  //add C ion
  if(ion_type == C_ION){
    //set mass
    mass = mass_matrix[cleavage_idx];

    if(constraint->mass_type == MONO){
      mass += MASS_NH3_MONO; 
    }
    else{ //average
      mass += MASS_NH3_AVERAGE; 
    }
    
    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, C_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge for C ion");
      return FALSE;
    }
  }
  
  //add X ion
  if(ion_type == X_ION){
    //set mass 
    mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[cleavage_idx];
    
    if(constraint->mass_type == MONO){
      mass += MASS_CO_MONO + MASS_H2O_MONO;     
    }
    else{ //average
      mass += MASS_CO_AVERAGE + MASS_H2O_AVERAGE; 
    }

    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, X_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge for X ion");
      return FALSE;
    }
  }

  //add Y ion
  if(ion_type == ALL_ION || ion_type == Y_ION){
    //set mass 
    mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[cleavage_idx];

    if(constraint->mass_type == MONO){
      mass += MASS_H2O_MONO; 
    }
    else{ //average
      mass += MASS_H2O_AVERAGE;
    }
    
    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, Y_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge Y ion");
      return FALSE;
    }
    
  }
  
  //add Z ion
  if(ion_type == Z_ION){
    //set mass 
    mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[cleavage_idx];

    if(constraint->mass_type == MONO){
      mass = mass - MASS_NH3_MONO + MASS_H2O_MONO; 
    }
    else{ //average
      mass = mass - MASS_NH3_AVERAGE + MASS_H2O_AVERAGE;
    }
    
    //add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, cleavage_idx, modification_counts, Z_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge Z ion");
      return FALSE;
    }
  }

  return TRUE;
}

/**
 * The engine of ion series, predicts all the ions from the parent peptide that
 * meet the ion constraint. All predicted ions are stored in the ion_series as ion objects.
 * Predict ion series
 */
void predict_ions(
  ION_SERIES_T* ion_series ///< the ion series to predict ions for
  )
{
  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  int modification_counts[MAX_MODIFICATIONS];
  int* modifications = constraint->modifications;
  
  //modification indecies
  int nh3_idx = 0;
  int h2o_idx = 0;
  int isotope_idx = 0;
  int flank_idx = 0;

  //determine the increment value for index (either 1 or -1)
  int nh3_increment = 1;
  int h2o_increment = 1;

  //reset modification increment
  if(modifications[NH3] < 0){
    nh3_increment = -1;
  }
  if(modifications[H2O] < 0){
    h2o_increment = -1;
  }


  //create a mass matrix
  float* mass_matrix = create_ion_mass_matrix(ion_series->peptide, constraint->mass_type);  
  
  //get peptide length
  int peptide_length = (int)mass_matrix[0];

  //add ions for each cleavage indecies
  int cleavage_idx = 1;
  for(; cleavage_idx < peptide_length; ++cleavage_idx){
    //FIXME add additional for loops when add new modiications
    //iterate over all different modification counts
    for(nh3_idx = 0; abs(nh3_idx) <= abs(modifications[NH3]); nh3_idx += nh3_increment){
      modification_counts[NH3] = nh3_idx;
      for(h2o_idx = 0; abs(h2o_idx) <= abs(modifications[H2O]); h2o_idx += h2o_increment){
        modification_counts[H2O] = h2o_idx;
        for(isotope_idx = 0; isotope_idx <= modifications[ISOTOPE]; ++isotope_idx){
          modification_counts[ISOTOPE] = isotope_idx;
          for(flank_idx = 0; flank_idx <= modifications[FLANK]; ++flank_idx){
            modification_counts[FLANK] = flank_idx;
            
            //add ions
            if(!add_ions(ion_series, cleavage_idx, mass_matrix, modification_counts)){
              carp(CARP_FATAL, "failed to add ion to ion series");
              free(ion_series);
              exit(1);
            }
          }
        }
      }
    }
  }
  //set predicted to TRUE
  ion_series->is_predicted = TRUE;
  
  //free mass matrix
  free(mass_matrix);
}


/**
 * Copies ion_series object from src to dest.
 *  must pass in a memory allocated ION_SERIES_T* dest
 */
void copy_ion_series(
  ION_SERIES_T* src,///< ion to copy from -in
  ION_SERIES_T* dest///< ion to copy to -out
  )
{
  ION_T* src_ion = NULL;
  ION_T* dest_ion = NULL;
  
  dest->peptide = my_copy_string(src->peptide);
  dest->charge = src->charge;

  //add copy of pointer ion constraint
  dest->constraint = src->constraint;

  //add copy ion, add ion_filtered_iterator
  ION_ITERATOR_T* iterator = new_ion_iterator(src);

  //iterate over all ions in src and copy them into dest
  while(ion_iterator_has_next(iterator)){
    src_ion = ion_iterator_next(iterator);
    //add ion
    dest_ion = allocate_ion(); 
    copy_ion(src_ion, dest_ion, dest->peptide);
    dest->ions[dest->num_ions++] = dest_ion;
  }
  //free up iterator
  free_ion_iterator(iterator);

  dest->is_predicted = TRUE;
}

/*************************************
 * ION_SERIES_T: get and set methods
 ************************************/

/**
 * User should not free the peptide sequence seperate from the ion_series
 *\returns a pointer to the original parent peptide sequence of the ion_series object
 */
char* get_ion_series_peptide(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  )
{
  return ion_series->peptide;
}

/**
 * copies in the peptide sequence to heap allocated sequence.
 * set the parent peptide sequence of the ion_series object
 */
void set_ion_series_peptide(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  char* peptide///< the peptide sequence to set -in
  )
{
  ion_series->peptide = my_copy_string(peptide);
}

/**
 *\returns the charge of the ion_series object
 */
int get_ion_series_charge(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  )
{
  return ion_series->charge;
}

/**
 * set the charge of the ion_series object
 */
void set_ion_series_charge(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  int charge///< the charge of the ion -in
  )
{
  ion_series->charge = charge;
}

/**
 *\returns the constraint of the ion_series object
 */
ION_CONSTRAINT_T* get_ion_series_ion_constraint(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  )
{
  return ion_series->constraint;
}

/**
 * frees the old ion_constraint, and replace with the given new constraint 
 * set the of the ion_series object
 */
void set_ion_series_ion_constraint(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  ION_CONSTRAINT_T* constraint///<  -in
  )
{
  free_ion_constraint(ion_series->constraint);
  ion_series->constraint = constraint;
}


/*************************
 * ION_CONSTRAINT methods
 *************************/

/**
 *\returns an empty heap allocated ion_constraint
 */
ION_CONSTRAINT_T* allocate_ion_constraint(void){
  ION_CONSTRAINT_T* constraint = (ION_CONSTRAINT_T*)mycalloc(1, sizeof(ION_CONSTRAINT_T));
  return constraint;
}

/**
 * modification, add more modifications as needed
 * copies the modifications, into its own array only if use_neutral_losses == TURE
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint(
  BOOLEAN_T use_neutral_losses, ///< A boolean to determine if the ions series should include neutral losses
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  int nh3_count, ///< the number of modifications of nh3
  int h2o_count, ///< the number of modifications of h2o
  int isotope_count, ///< the number of modifications of isotope
  int flank_count ///< the number of modifications of flank
  )
{
  ION_CONSTRAINT_T* constraint = allocate_ion_constraint();
  constraint->use_neutral_losses = use_neutral_losses;

  //if use natural loss, copy
  if(use_neutral_losses){
    //set all modifications count
    constraint->modifications[NH3] = nh3_count;
    constraint->modifications[H2O] = h2o_count;
    constraint->modifications[ISOTOPE] = isotope_count;
    constraint->modifications[FLANK] = flank_count;
  }
  
  constraint->mass_type = mass_type;
  constraint->max_charge = max_charge;
  constraint->ion_type = ion_type;

  return constraint;
}

/**
 * Frees an allocated ion_constraint object.
 */
void free_ion_constraint(
  ION_CONSTRAINT_T* ion_constraint///< the ion constraints to enforce -in
  )
{
  free(ion_constraint);
}

/**
 * copies ion_constraint object from src to dest
 * must pass in a memory allocated ION_CONSTRAINT_T dest
 */
void copy_ion_constraint(
  ION_CONSTRAINT_T* src,///< ion_constraint to copy from -in
  ION_CONSTRAINT_T* dest///< ion_constraint to copy to -out
  )
{
  int modification_idx = 0;
  dest->use_neutral_losses = src->use_neutral_losses;

  //if use natural loss, copy
  if(src->use_neutral_losses){
    //iterate over all modifications a update new constraint
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      dest->modifications[modification_idx] = src->modifications[modification_idx];
    }
  }
  
  dest->mass_type = src->mass_type;
  dest->max_charge = src->max_charge;
  dest->ion_type = src->ion_type;
}

/** 
 * Determines if a ion satisfies a ion_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T ion_constraint_is_satisfied(
   ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
   ION_T* ion ///< query ion -in
   )
{
  int* counts = NULL;
  int modification_idx = 0;

  //check ion type
  if(get_ion_type(ion) != ion_constraint->ion_type &&
     !((ion_constraint->ion_type == ALL_ION) && 
       (get_ion_type(ion) == B_ION || get_ion_type(ion) == Y_ION))){
       
    return FALSE;
  }
  
  //check charge
  if(get_ion_charge(ion) > ion_constraint->max_charge){
    return FALSE;
  }
  
  //check modifications
  counts = get_ion_modification_counts(ion);
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    if(counts[modification_idx] != ion_constraint->modifications[modification_idx]){
      return FALSE;
    }
  }
  
  //FIXME, add more checks here as more contraints are added

  return TRUE;
}


/**************************
 *  ION_ITERATOR_T object
 **************************/

/**
 * Instantiates a new ion_iterator object from ion_series.
 * \returns a ION_ITERATOR_T object.
 */
ION_ITERATOR_T* new_ion_iterator(
  ION_SERIES_T* ion_series ///< ion_series to iterate -in
  )
{
  ION_ITERATOR_T* iterator = (ION_ITERATOR_T*)mycalloc(1, sizeof(ION_ITERATOR_T));
  iterator->ion_series = ion_series;
  return iterator;
}

/**
 * does not free ions
 * Frees an allocated ion_iterator object.
 */
void free_ion_iterator(
  ION_ITERATOR_T* ion_iterator///< free ion_iterator -in
  )
{
  free(ion_iterator);
}

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T ion_iterator_has_next(
  ION_ITERATOR_T* ion_iterator///< is there a next ion? -in
  )
{
  return ion_iterator->ion_series->num_ions < ion_iterator->ion_idx;
}

/**
 * The basic iterator function next.
 */
ION_T* ion_iterator_next(
  ION_ITERATOR_T* ion_iterator///< return the next ion -in
  )
{
  ++ion_iterator->ion_idx;
  
  return ion_iterator->ion_series->ions[ion_iterator->ion_idx-1]; 
}

/**********************************
 * ION_FILTERED_ITERATOR_T object
 **********************************/

/**
 * sets up the iterator for next iteration.
 * 
 *\returns TRUE if successfully sets up the ion_filtered_iterator for next iteration
 */
BOOLEAN_T setup_ion_filtered_iterator(
  ION_FILTERED_ITERATOR_T* ion_iterator///< free ion_iterator -in
  )
{
  BOOLEAN_T found = FALSE;
  ION_T* ion = NULL;

  //iterate over ions until discovers the first ion that meets the ion constraint
  while(ion_iterator->ion_idx < ion_iterator->ion_series->num_ions){
    //get next ion
    ion = ion_iterator->ion_series->ions[ion_iterator->ion_idx];
    
    //check if the current ion satisfies the ion_constraint for the iterator
    if(ion_constraint_is_satisfied(ion_iterator->constraint, ion)){
      found = TRUE;
      ion_iterator->ion = ion;
      ++ion_iterator->ion_idx;
      break;
    }
    ++ion_iterator->ion_idx;
  }
  
  ion_iterator->has_next = found;

  return TRUE;
}

/**
 * Only copies in the constraint as pointer
 * Instantiates a new ion_filtered_iterator object from ion_series.
 * \returns a ION_FILTERED_ITERATOR_T object.
 */
ION_FILTERED_ITERATOR_T* new_ion_filtered_iterator(
  ION_SERIES_T* ion_series, ///< ion_series to iterate -in
  ION_CONSTRAINT_T* constraint  ///< ion_constraint which returned ions satisfy
  )
{
  ION_FILTERED_ITERATOR_T* iterator = 
    (ION_FILTERED_ITERATOR_T*)mycalloc(1, sizeof(ION_FILTERED_ITERATOR_T));
  
  //set constraint, ion_series
  iterator->constraint = constraint;
  iterator->ion_series = ion_series;
  iterator->has_next = FALSE;

  //initialize iterator
  setup_ion_filtered_iterator(iterator);

  return iterator;
}        

/**
 * The constraint is NOT freed from the iterator.
 * Frees an allocated ion_filtered_iterator object.
 */
void free_ion_filtered_iterator(
  ION_FILTERED_ITERATOR_T* ion_iterator///< free ion_iterator -in
  )
{
  free(ion_iterator);
}

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T ion_filtered_iterator_has_next(
  ION_FILTERED_ITERATOR_T* ion_iterator///< is there a next ion? -in
  )
{
  return ion_iterator->has_next;
}

/**
 * The basic iterator function next.
 */
ION_T* ion_filtered_iterator_next(
  ION_FILTERED_ITERATOR_T* ion_iterator///< return the next ion -in
  )
{
  ION_T* next_ion = NULL;
  
  //check if a ion is present to return
  if(ion_iterator->has_next){
    carp(CARP_FATAL, "index out of bounds for ion_filtered_iterator");
    free(ion_iterator);
    exit(1);
  }
  
  next_ion = ion_iterator->ion;
  
  //re-initialize iterator
  setup_ion_filtered_iterator(ion_iterator);

  return next_ion;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
