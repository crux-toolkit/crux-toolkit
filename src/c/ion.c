/*****************************************************************************
 * \file ion.c
 * $Revision: 1.8 $
 * \brief: Object for representing a single ion.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "objects.h"
#include "ion.h"
#include "mass.h"
#include "utils.h"

// MAX_MODIFICATIONS 4, defined in ion.h, because ion_series.c needs the information

/**
 * Array to store the modification masses
 */
float modification_masses[MAX_MODIFICATIONS];

/**
 * Have we initialized the modification_masses?
 */
BOOLEAN_T initialized_modification_masses = FALSE;

/**
 * \struct ion
 * \brief An object for representing a (fragment) ion of a peptide.
 */
struct ion {
  ION_TYPE_T type;  ///< type of the ion 
  int cleavage_idx; ///< index of peptide amide that fragments to form this ion, starting from the N-term end 
  // N.b. this is different than the b1,y1 index, in that it always starts
  // from the N-term
  int charge; ///< the ion charge
  char* peptide_sequence; ///< the peptide sequence that fragments to form this ion
  int modification_counts[MAX_MODIFICATIONS]; ///< an array of the number of different ion modifications
  float ion_mass_z;   ///< The mass/z of the ion. 
};


/**
 * initializes the mass array
 */
void initialize_modification_masses(
   MASS_TYPE_T mass_type ///< mass type (average, mono) -in
)
{
  //set modification mass
  if(mass_type == MONO){
    modification_masses[NH3] = MASS_NH3_MONO;
    modification_masses[H2O] = MASS_H2O_MONO ;
    modification_masses[ISOTOPE] = 1; //FIXME check this!!!
    modification_masses[FLANK] = 1; //FIXME check this!!!
  }
  else if(mass_type == AVERAGE){
    modification_masses[NH3] = MASS_NH3_AVERAGE;
    modification_masses[H2O] = MASS_H2O_AVERAGE;
    modification_masses[ISOTOPE] = 1; //FIXME check this!!!
    modification_masses[FLANK] = 1; //FIXME check this!!!
  }

  initialized_modification_masses = TRUE;
}


/**
 * \returns An (empty) ion object.
 */
ION_T* allocate_ion(void){
  ION_T* ion = (ION_T*)mycalloc(1, sizeof(ION_T));
  return ion;
}

/**
 * helper function
 * only copies the pointer to the peptide sequence
 * creates an heap allocated ion, ion mass is not calculated.
 * \returns an ION_T object
 */
ION_T* new_basic_ion (
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide ///< location for the new ion -in
  )
{
  //allocate new ion
  ION_T* ion = (ION_T*)mycalloc(1, sizeof(ION_T));
  ion->type = type;
  ion->cleavage_idx = cleavage_idx;
  ion->charge = charge;
  ion->peptide_sequence = peptide;
  return ion;
}


/**
 * only copies the pointer to the peptide sequence
 * ion mass/z with out any modification
 * \returns an ION_T object
 */
ION_T* new_ion (
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  //get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  //calculate and set ion mass/z
  if(!calc_ion_mass_z(ion, mass_type, FALSE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
    exit(1);
  }
  
  return ion;
}

/**
 * only copies the pointer to the peptide sequence
 * inputs a array of all the modification counts
 * \returns an ION_T object
 */
ION_T* new_modified_ion(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  int* modification_counts ///< an array of modification counts for each modification -in
  )
{
  //get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  //set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = modification_counts[modification_idx];
  }
  
  //calculate and set ion mass/z
  if(!calc_ion_mass_z(ion, mass_type, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
    exit(1);
  }
  
  return ion;
}

/**
 * only copies the pointer to the peptide sequence
 * inputs a array of all the modification counts
 * inputs the pre modified mass, of just all AA mass summed up.
 * \returns an ION_T object
 */
ION_T* new_modified_ion_with_mass(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  float base_mass, ///< the base mass of the ion -in
  int* modification_counts ///< an array of modification counts for each modification -in
  )
{
  //get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  //set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = modification_counts[modification_idx];
  }
  
  //calculate and set ion mass/z
  if(!calc_ion_mass_z_with_mass(ion, mass_type, base_mass, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
    exit(1);
  } 
  return ion;
}

/**
 * only copies the pointer to the peptide sequence
 * inputs the pre modified mass, of just all AA mass summed up.
 * \returns an ION_T object
 */
ION_T* new_ion_with_mass(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  float base_mass ///< the base mass of the ion -in
  )
{
  //get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  //set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = 0;
  }
  
  //calculate and set ion mass/z
  if(!calc_ion_mass_z_with_mass(ion, mass_type, base_mass, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
    exit(1);
  } 
  return ion;
}


/**
 * frees A ION_T object
 */
void free_ion (
  ION_T* garbage_ion ///< the ion to free -in
  )
{
  free(garbage_ion);
}

/**
 * prints the location and fields of ION_T object to the file, in the
 * following format:
 *
 * m/z \\t mass \\t charge \\t ion-series \\t  ...
 *  peptide-bond-index \\t modifications \n
 *
 * Where:
 *
 * m/z - is the ion's mass-to-charge
 * mass - is the ion's (charged) mass
 * charge - is the ion's charge e.g. 1,2,3
 * ion-series - is one of (b,y,p)
 * bond-index - is in [1...n), where n is peptide length
 * modifications - is one of (none|nh3|h2o)
 *
 * if the ion has more than one modification, each will be printed on a
 * separate line, with the necessary number of tabs to right justify
 */
void print_ion(
  ION_T* ion, ///< print this ion -in
  FILE* file ///< to this file -in
  )
{
  //print all fields of ion
  fprintf(file, "%.2f\t%.2f\t%d\t%d\t%d", ion->ion_mass_z, (ion->ion_mass_z)*ion->charge, ion->charge, 
          (int)ion->type, ion->cleavage_idx);

  //iterate over all modification counts
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    fprintf(file,"\t%d", ion->modification_counts[modification_idx]);
  }
  
  fprintf(file,"\n");
}

/**
 *\return the modified mass_z accodring to the modification type
 */
float modify_ion_mass_z(
  float mass_z, ///< the pre-modified ion mass_z -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  int charge, ///< charge of the ion
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  //initalize the modified masses(average|mono);
  if(!initialized_modification_masses){
    initialize_modification_masses(mass_type);
  }
  
  return  mass_z + (modification_masses[(int)ion_modification]/(float)charge) * modification_count;  
}


/**
 * Adds the given ION_MODIFICATION to this ion
 */
void add_modification(
  ION_T* ion,///< ion to which to add the modification -mod
  ION_MODIFICATION_T modification, ///< add this modification to the ion -in
  int modification_count, ///< the number of modifications -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  //update modification count
  ion->modification_counts[(int)modification] += modification_count;  
  //reset ion mass_z
  ion->ion_mass_z =  modify_ion_mass_z(ion->ion_mass_z, modification_count, modification, 
                                       ion->charge, mass_type);
}


/**
 *\returns the ion's AA mass added all up
 */
float get_ion_mass(
  ION_T* ion, ///< the working ion -out/in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )                   
{
  char* ion_sequence = NULL;
  ION_TYPE_T ion_type = ion->type;
  int ion_length = 0;
  BOOLEAN_T memory_used = FALSE;
  float mass = 0;
  BOOLEAN_T reverse = FALSE;

  //get sequence for x,y,z ion
  if(ion_type == X_ION ||ion_type == Y_ION || ion_type == Z_ION){
    //convert the cleavage index into the actually index that start from the left.
    int real_cleavage_idx = strlen(ion->peptide_sequence) - ion->cleavage_idx;
    ion_sequence = &(ion->peptide_sequence[real_cleavage_idx]);
    ion_length = strlen(ion_sequence);
    reverse = TRUE;
  }
  //get sequence for a,b,c ion
  else{
    ion_length = ion->cleavage_idx;
    ion_sequence = mycalloc(ion_length+1, sizeof(char));
    strncpy(ion_sequence, ion->peptide_sequence, ion_length);
    memory_used = TRUE;
  }

  //add up all AA mass
  int ion_idx = 0;
  for(; ion_idx < ion_length; ++ion_idx){
    mass += get_mass_amino_acid(ion_sequence[ion_idx], mass_type);
  }
  
  //free ion sequence, only if memory allocated (a,b,c ions);
  if(memory_used){
    free(ion_sequence);
  }

  //if X,Y,Z ion add H2O
  if(reverse){
    if(mass_type == AVERAGE){
      return mass + MASS_H2O_AVERAGE;
    }
    return mass + MASS_H2O_MONO;
  }
  
  return mass;
}

/**
 *\return the modified mass accodring to the modification type
 */
float modify_ion_mass(
  float mass, ///< the pre-modified ion mass -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  //initalize the modified masses(average|mono);
  if(!initialized_modification_masses){
    initialize_modification_masses(mass_type);
  }
  
  return  mass + modification_masses[(int)ion_modification] * modification_count;
}

/**
 * is_modified, indiciates if there are any modification to the ion
 * speeds up the proccess if FALSE.
 *\returns TRUE if successfully computes the mass/z of the ion, else FALSE
 */
BOOLEAN_T calc_ion_mass_z_with_mass(
  ION_T* ion, ///< the working ion -out/in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  float mass, ///< the basic mass of the ion -in
  BOOLEAN_T is_modified ///< are there any modifications for this ion? -in
  )
{
  float h_mass = MASS_H_MONO;

  //alter mass according to the modification
  if(is_modified){
    //iterate over all type of modifications for ion
    int modification_idx = 0;
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      //update ion mass if modification is needed
      if(ion->modification_counts[modification_idx] != 0){
        mass = modify_ion_mass(mass, ion->modification_counts[modification_idx], 
                               (ION_MODIFICATION_T)modification_idx, mass_type);
      }
    }
  }
  
  //reset h mass if needed
  if(mass_type == AVERAGE){
    h_mass = MASS_H_AVERAGE;
  }

  //convert mass to m/z and assigned to ion
  ion->ion_mass_z = (mass + (h_mass*(float)ion->charge))/(float)ion->charge;

  return TRUE;
}

/**
 * is_modified, indiciates if there are any modification to the ion
 * speeds up the proccess if FLASE.
 *\returns TRUE if successfully computes the mass/z of the ion, else FALSE
 */
BOOLEAN_T calc_ion_mass_z(
  ION_T* ion, ///< the working ion -out/in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  BOOLEAN_T is_modified ///< are there any modifications for this ion? -in
  )
{
  float mass = get_ion_mass(ion, mass_type);
  return calc_ion_mass_z_with_mass(ion, mass_type, mass, is_modified);
}


/**
 * Copies ion object from src to dest.
 * must pass in a seperate pointer peptide sequence from its own ion_series object
 * must pass in a memory allocated ION_T* dest
 */
void copy_ion(
  ION_T* src,///< ion to copy from -in
  ION_T* dest,///< ion to copy to -out
  char* peptide_sequence ///< the peptide sequence that the dest should refer to -in
  )
{
  dest->type = src->type;
  dest->cleavage_idx = src->cleavage_idx;
  dest->charge = src->charge;
  
  //set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    dest->modification_counts[modification_idx] = src->modification_counts[modification_idx];
  }
  
  dest->ion_mass_z = src->ion_mass_z;
  dest->peptide_sequence = peptide_sequence;
}


/**
 *\returns TRUE if forward ion_type(A,B,C), else reverse ion_type(X,Y,Z) FALSE
 */
BOOLEAN_T is_forward_ion_type(
  ION_T* ion ///< the ion to check if can lose nh3 -in                         
  )
{
  //is ion forward type?
  if(ion->type == B_ION ||
     ion->type == A_ION ||
     ion->type == C_ION)
    {
      return TRUE;
    }

  //reverse type ion
  return FALSE;
}

/**
 *\returns TRUE if the ion has modifications, else FALSE
 */
BOOLEAN_T ion_is_modified(
  ION_T* ion ///< the ion to check if can lose nh3 -in
  )
{
  int by_modification = 0;

  //only add ions with no modifications
  for(; by_modification < MAX_MODIFICATIONS; ++by_modification){
    if(ion->modification_counts[by_modification] != 0){
      return TRUE;
    }
  }
  
  return FALSE;
}

/*********************************
 * get, set methods for ion fields
 *********************************/

/**
 * \returns the location of ION_T object
 */
float get_ion_mass_z(
  ION_T* working_ion///< return the location of this ion -in 
  )
{
  return working_ion->ion_mass_z;
}


/**
 * sets the mass/z of the ION_T object
 * while this can be calculated from the char*, cleavage_idx and
 * modifications, it allows some optimizations if we allow it to be set
 * instead
 */
void set_ion_mass_z(
  ION_T* working_ion, ///<set the m/z location of this ion -out
  float mass_z ///< the m/z location -in
  )
{
  working_ion->ion_mass_z = mass_z;
}

/**
 * return the cleavage_idx of the ion object
 */
int get_ion_cleavage_idx(
  ION_T* working_ion///< the working ion -in                          
  )
{
  return working_ion->cleavage_idx;
}

/**
 * set the cleavage_idx of the ion object
 */
void set_ion_cleavage_idx(
  ION_T* working_ion, ///< the working ion -out
  int cleavage_idx ///< the cleavage index in the peptide -in
  )
{
  working_ion->cleavage_idx = cleavage_idx;
}

/**
 * return the charge of the ion object
 */
int get_ion_charge(
  ION_T* working_ion ///< the working ion -in                          
  )
{
  return working_ion->charge;
}

/**
 * set the charge of the ion object
 */
void set_ion_charge(
  ION_T* working_ion, ///< the working ion -out
  int charge ///< the charge of this ion -in
  )
{
  working_ion->charge = charge;
}

/**
 * return the ION_TYPE_T of the ion object
 */
ION_TYPE_T get_ion_type(
  ION_T* working_ion ///< the working ion -in                          
  )
{
  return working_ion->type;
}

/**
 * set the ION_TYPE_T of the ion object
 */
void set_ion_type(
  ION_T* working_ion, ///< the working ion -out
  ION_TYPE_T ion_type ///< the ion type of this ion -in 
  )
{
  working_ion->type = ion_type;
}

/**
 * return the parent peptide sequence of the ion object
 * returns a pointer to the sequence, should not free
 */
char* get_ion_peptide_sequence(
  ION_T* working_ion ///< the working ion -in                          
  )
{
  return working_ion->peptide_sequence;
}

/**
 * return a pointer to the modification_count array of the ion object
 */
int* get_ion_modification_counts(
  ION_T* working_ion ///< the working ion -in                          
  )
{
  return working_ion->modification_counts;
}

/**
 * return the count of in the modification_count array of the ion object
 */
int get_ion_single_modification_count(
  ION_T* working_ion, ///< the working ion -in                          
  ION_MODIFICATION_T mod_type ///< the modification count wanted -in
  )
{
  return working_ion->modification_counts[mod_type];
}


/**
 * set the parent peptide_sequence of the ion object
 */
void set_ion_peptide_sequence(
  ION_T* working_ion, ///< the working ion -out
  char* peptide_sequence ///< the parent peptide's sequence of this ion -in 
  )
{
  working_ion->peptide_sequence = peptide_sequence;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

