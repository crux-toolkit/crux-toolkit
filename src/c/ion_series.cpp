/*************************************************************************//**
 * \file ion_series.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 21 Sep 2006
 * DESCRIPTION: code to support working with a series of ions
 * REVISION: $Revision: 1.52 $
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
#include "parameter.h"
#include "peptide.h"
#include "mass.h"
#include "spectrum.h"

#define BINARY_GMTK 1
#define PRINT_NULL_IONS 1
#define MIN_FRAMES 3
#define MAX_IONS 10000
#define MAX_NUM_ION_TYPE 8 // number of different ion_types

/**
 * \struct ion_series
 * \brief An object to represent a series of ions, and organize them.
 * For which additional data structures will be created as needed 
 * loss_limit can be equal to NULL, thus if need to use should always
 * check that it is not NULL.
 */
struct ion_series {
  // TODO change name to unmodified_char_seq
  char* peptide; ///< The peptide sequence for this ion series
  MODIFIED_AA_T* modified_aa_seq; ///< sequence of the peptide
  FLOAT_T peptide_mass; ///< The peptide neutral mass. For efficiency. 
  int charge; ///< /<The charge state of the peptide for this ion series
  ION_CONSTRAINT_T* constraint; ///< The constraints which these ions obey
  ION_T* ions[MAX_IONS]; ///< The ions in this series
  int num_ions; ///< the number of ions in this series
  BOOLEAN_T is_predicted; ///< has this ion_series been predicted already?
  int num_specific_ions[MAX_NUM_ION_TYPE]; 
    ///< the number of ions of a specific ion_type
  ION_T* specific_ions[MAX_NUM_ION_TYPE][MAX_IONS]; 
    ///< specific ions in the series, reference to master array of ions
  LOSS_LIMIT_T* loss_limit; 
    ///< nh3, h2o loss limit for a given cleavage index, 
    ///< before using this array should always sheck if not NULL
  int peptide_length;   ///< the length of the peptide
};
// ??? what is the difference between peptide_length and num_ions

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
 * \struct ion_constraint
 * \brief An object to represent the contraints which the ions in this
 * series obey.
 */
struct ion_constraint {
  BOOLEAN_T use_neutral_losses; ///< Should ions include neutral losses
  int modifications[MAX_MODIFICATIONS]; 
    ///< an array to indicate which modifications to perform
  MASS_TYPE_T mass_type; 
    ///< the mass_type to use MONO|AVERAGE
  int max_charge; 
    ///< maximum charge of the ions, cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type; 
    ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion; 
    ///< does a precursor-ion satisfy this constraint
  int min_charge; 
  BOOLEAN_T exact_modifications; 
    ///< TRUE  = ints in modfications array indicate exact number of mods
    ///< FALSE = ints in modfications array indicate maximum number of mods
  unsigned int pointer_count; ///< Number of pointers referencing me 
};

/**
 * \struct ion_iterator
 * \brief An object to iterate over all ion objects in the ion_series
 */
struct ion_iterator {
  ION_SERIES_T* ion_series; ///< the ion series that the ion we are iterating
  int ion_idx; ///< the current ion that is being returned 
};

/**
 * \struct ion_filtered_iterator
 * \brief An object to iterate over ion objects that meet constraint in
 * the ion_series 
 */
struct ion_filtered_iterator {
  ION_SERIES_T* ion_series; ///< the ion series that the ion we are iterating
  ION_CONSTRAINT_T* constraint; ///< constraints which the ions obey
  BOOLEAN_T has_next; ///< the boolean which the iterator has a next ion
  int ion_idx; ///< the current ion that is being returned 
  ION_T* ion; ///< the next ion to return when called upon
  ION_T** ion_array; 
  ///< the specfic ion array we are iterating over, B ion, Y ion or all
  int array_size; ///< size of the ion_array
};

/**
 * \returns An (empty) ion_series object.
 */
ION_SERIES_T* allocate_ion_series(void){
  int ion_type_idx = 0;

  ION_SERIES_T* ion_series = (ION_SERIES_T*)mycalloc(1, sizeof(ION_SERIES_T));
  ion_series->is_predicted = FALSE;
  
  // initialize all num_specific_ion count to 0
  for(ion_type_idx = 0; ion_type_idx < MAX_NUM_ION_TYPE; ++ion_type_idx){
    ion_series->num_specific_ions[ion_type_idx] = 0;
  }

  return ion_series;
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
ION_SERIES_T* new_ion_series(
  const char* peptide, ///< The peptide for this ion series. -in
  int charge, ///< The charge for this ion series -in
  ION_CONSTRAINT_T* constraint ///< constraints which these ions obey.
  )
{
  ION_SERIES_T* ion_series = allocate_ion_series();

  // copy the peptide sequence
  ion_series->peptide = my_copy_string(peptide);
  ion_series->modified_aa_seq = convert_to_mod_aa_seq(peptide);
  ion_series->peptide_mass = calc_sequence_mass(peptide, MONO);
  ion_series->charge = charge;
  ion_series->constraint = constraint;
  ion_series->peptide_length = strlen(peptide);
  
  // create the loss limit array
  ion_series->loss_limit = 
    (LOSS_LIMIT_T*)mycalloc(ion_series->peptide_length, sizeof(LOSS_LIMIT_T));

  return ion_series;
}

/**
 * \brief Creates an ion_series object without a specific peptide.
 * Peptide details are added with the "update_ion_series" method.
 * \returns A newly allocated ion_series object that must be updated
 * for each peptide instance.
 */
ION_SERIES_T* new_ion_series_generic(
  ION_CONSTRAINT_T* constraint, ///< constraints which these ions obey.
  int charge ///< The charge for this ion series -in
  )
{
  ION_SERIES_T* ion_series = allocate_ion_series();
  ion_series->constraint = constraint;
  ion_series->charge = charge;
  // use max peptide len so loss_limit array can be used for any peptide
  ion_series->loss_limit = 
    (LOSS_LIMIT_T*)mycalloc(get_int_parameter("max-length"), 
                            sizeof(LOSS_LIMIT_T));
  return ion_series;
}

/**
 * \brief Updates an ion_series to a specific instance of a peptide
 * sequence. If the ion_series has already generated ions, they will
 * be free-ed. A copy of the peptide sequence is made and all other
 * variables (except ion_constraint) are also updated for the new
 * peptide sequence. 
 */
void update_ion_series(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  char* peptide, ///< The peptide sequence with no mod characters. -in
  MODIFIED_AA_T* mod_seq ///< modified version of char* sequence -in
  ) 
{
  int ion_type_idx = 0;

  // Initialize the ion_series object for the new peptide sequence
  
  // free old peptide sequence
  if(ion_series->peptide){
    free(ion_series->peptide);
  }
  if(ion_series->modified_aa_seq){
    free(ion_series->modified_aa_seq);
  }
  
  // iterate over all ions, and free them
  while(ion_series->num_ions > 0){
    free_ion(ion_series->ions[ion_series->num_ions-1]);
    --ion_series->num_ions;
  }
  
  // initialize all num_specific_ion count back to 0
  int ion_idx;
  for(ion_idx=0; ion_type_idx < MAX_NUM_ION_TYPE; ++ion_type_idx){
    ion_series->num_specific_ions[ion_type_idx] = 0;
  }
  
  ion_series->num_ions = 0;
  ion_series->is_predicted = FALSE;
  
  // set ion_series for new instance of peptide
  
  // copy the peptide sequence
  ion_series->peptide = my_copy_string(peptide);
  ion_series->peptide_length = strlen(peptide);
  ion_series->modified_aa_seq = copy_mod_aa_seq(mod_seq, 
                                                ion_series->peptide_length);
  
  // Initialize the loss limit array for the new peptide
  for(ion_idx =0; ion_idx < ion_series->peptide_length; ++ion_idx){
    ion_series->loss_limit->h2o = 0;
    ion_series->loss_limit->nh3 = 0;
    // add more initialize count if more added
  }
}


/**
 * Frees an allocated ion_series object.
 */
void free_ion_series(
  ION_SERIES_T* ion_series ///< the ion collection to free - in
  )
{
  if(ion_series->peptide){
    free(ion_series->peptide);
  }
  if(ion_series->modified_aa_seq){
    free(ion_series->modified_aa_seq);
  }
  if(ion_series->loss_limit){
    free(ion_series->loss_limit);
  }
  // free constraint?

  // iterate over all ions, and free them
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
  // check if the ions has been already predicted
  if(!ion_series->is_predicted){
    carp(CARP_ERROR, "ion series has not predicted ions");
    return;
  }
  
  // print header
  fprintf(file, "m/z\tmass\tcharge\tion-series\tpeptide-bond-index\tNH3\tH2O\tISOTOPE\tFLANK\n");
  
  
  // print each ion in the ion series
  int ion_idx;
  for(ion_idx=0; ion_idx < ion_series->num_ions; ++ion_idx){
    print_ion(ion_series->ions[ion_idx], file);
  }
}

/**
 * Prints a ion_series object to file, in GMTK single-ion format.
 */
void print_ion_series_single_gmtk(
  ION_SERIES_T* ion_series,         ///< ion_series to print -in 
  ION_CONSTRAINT_T* ion_constraint, ///< ion_constraint to obey -in 
  FILE* file,                       ///< file output
  int sentence_idx){

  // create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = 
    new_ion_filtered_iterator(ion_series, ion_constraint);
  
  // foreach ion in ion iterator, add matched observed peak intensity
  ION_T* ion;
  int frame_idx = 0;
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ion = ion_filtered_iterator_next(ion_iterator);
    
#ifdef BINARY_GMTK
    print_ion_gmtk_single_binary(ion, file, sentence_idx, frame_idx);
#else
    print_ion_gmtk_single(ion, file);
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
  
  free_ion_filtered_iterator(ion_iterator);
}

/**
 * Prints a ion_series object to file, in GMTK paired-ion format.
 */
void print_ion_series_paired_gmtk(
  ION_SERIES_T* ion_series, ///< ion_series to print -in 
  ION_CONSTRAINT_T* first_ion_constraint, ///< ion_constraint to obey -in 
  ION_CONSTRAINT_T* second_ion_constraint, ///< ion_constraint to obey -in 
  FILE* file, ///< file output
  int sentence_idx
  )
{
  
  // create the filtered iterator that will select among the ions
  ION_FILTERED_ITERATOR_T* ion_iterator = 
    new_ion_filtered_iterator(ion_series, first_ion_constraint);
  
  // foreach ion in ion iterator, add matched observed peak intensity
  int frame_idx = 0;
  while(ion_filtered_iterator_has_next(ion_iterator)){
    ION_T* first_ion = ion_filtered_iterator_next(ion_iterator);
    int cleavage_idx = get_ion_cleavage_idx(first_ion);
    ION_T* second_ion = get_ion_series_ion(
                            ion_series, second_ion_constraint, cleavage_idx);
    if ( (first_ion == NULL) || (second_ion == NULL) ){
      continue;
    }
    print_ion_gmtk_paired_binary(
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
 
  free_ion_filtered_iterator(ion_iterator);
}


/**
 * \brief Find instances of amino acid which can incur neutral
 * losses: H2O (S|T|E|D), NH3(R|K|Q|N).  
 * Set the count of those observed so far for each cleavage index.
 * If no instance of amino acid, the count is assigned to 0
 * The information is used to determine if how many nh3 or h2o neutral
 * losses are possible. 
 */
void scan_for_aa_for_neutral_loss(
  ION_SERIES_T* ion_series ///< ion_series to print -in/out
  )
{
  int peptide_length = ion_series->peptide_length;
  char* sequence = ion_series->peptide;

  int h2o_aa = 0;
  int nh3_aa = 0;
  LOSS_LIMIT_T* loss_limit_count = NULL; // debug

  // search for the first instance of the amino acids
  int cleavage_idx;
  for(cleavage_idx=0; cleavage_idx < peptide_length; ++cleavage_idx){
    // is the AA  (S|T|E|D) ?
    if(sequence[cleavage_idx] == 'S' ||
       sequence[cleavage_idx] == 'T' ||
       sequence[cleavage_idx] == 'E' ||
       sequence[cleavage_idx] == 'D' )
      {
        loss_limit_count = &ion_series->loss_limit[cleavage_idx];
        loss_limit_count->h2o = ++h2o_aa;
        loss_limit_count->nh3 = nh3_aa;
      }
    // is the AA  (R|K|Q|N) ?
    else if(sequence[cleavage_idx] == 'R' ||
            sequence[cleavage_idx] == 'K' ||
            sequence[cleavage_idx] == 'Q' ||
            sequence[cleavage_idx] == 'N' )
      {
        loss_limit_count = &ion_series->loss_limit[cleavage_idx];
        loss_limit_count->nh3 = ++nh3_aa;
        loss_limit_count->h2o = h2o_aa;
      }
    else{
      loss_limit_count = &ion_series->loss_limit[cleavage_idx];
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
FLOAT_T* create_ion_mass_matrix(
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
void add_ion_to_ion_series(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -out
  ION_T* ion ///< ion to add -in
  )
{
  // add ion to ion series
  ion_series->ions[ion_series->num_ions++] = ion;   
  
  // add a pointer of ion to the specific ion_type array
  ion_series->specific_ions[get_ion_type(ion)]
    [ion_series->num_specific_ions[get_ion_type(ion)]++] = ion;
}

/**
 * helper function: add_ions
 * add all the ions to ion_series up to the max charge
 *\returns TRUE if successfully adds all ions, else FALSE
 */
BOOLEAN_T add_ions_by_charge(
  ION_SERIES_T* ion_series, ///< the ion series to predict ions for -in
  FLOAT_T mass, ///< the base mass of the ion to add
  int cleavage_idx, ///< the absolute cleavage index (A,B,C from left X,Y,Z from right)
  ION_TYPE_T ion_type ///< the ion type of the ions to be added
  )
{
  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  int charge_idx = 1;
  ION_T* ion = NULL;
  int max_charge;

  // check if there's enough space to add the new ions
  if(ion_series->num_ions + constraint->max_charge > MAX_IONS){
    carp(CARP_ERROR, "exceeds ion array size in ion_series");
    return FALSE;
  }
  
  // set the max charge, the maximum cannot exceed the precursor ion's charge
  if(constraint->max_charge > ion_series->charge){
    max_charge = ion_series->charge;
  }
  else{
    max_charge = constraint->max_charge;
  }

  // iterate over all different charge
  for(; charge_idx <= max_charge; ++charge_idx){
    // create ion
    ion = new_ion_with_mass(ion_type, 
                            cleavage_idx, 
                            charge_idx, 
                            ion_series->peptide, 
                            constraint->mass_type, 
                            mass); 
    // add ion to ion series
    add_ion_to_ion_series(ion_series, ion);
  }
  
  return TRUE;
}


/**
 * Creates all the ions with no modifications up to the max charge
 * Adds each ion to ion_series
 *\returns TRUE if successfully generates all the ions, else FALSE
 */
BOOLEAN_T generate_ions_no_modification(
  ION_SERIES_T* ion_series, ///< ion_series to modify -in/out
  FLOAT_T* mass_matrix ///< the mass matrix that stores the mass
  )
{
  if( ion_series == NULL || mass_matrix == NULL ){
    carp(CARP_ERROR,
         "Cannot generate ions from NULL ion series or mass matrix");
    return FALSE;
  }
  int cleavage_idx = 1;
  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  FLOAT_T mass = 0;

  // get peptide length
  int peptide_length = (int)mass_matrix[0];

  // iterate over all cleavage index
  for(; cleavage_idx < peptide_length; ++cleavage_idx){
    
    // add A ion
    if(constraint->ion_type == A_ION 
       || constraint->ion_type == BYA_ION 
       || constraint->ion_type == ALL_ION){

      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->mass_type == MONO){
        mass -= MASS_CO_MONO; 
      }
      else{ // average
        mass -= MASS_CO_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, A_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for A ion");
      return FALSE;
      }
    }
    
    // add B ion
    if(constraint->ion_type == ALL_ION 
       || constraint->ion_type == BY_ION
       || constraint->ion_type == BYA_ION
       || constraint->ion_type == B_ION){
      
      // set mass
      mass = mass_matrix[cleavage_idx];
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, B_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for B ion");
        return FALSE;
      }
    }
    
    // add C ion
    if(constraint->ion_type == C_ION || constraint->ion_type == ALL_ION){
      // set mass
      mass = mass_matrix[cleavage_idx];
      
      if(constraint->mass_type == MONO){
        mass += MASS_NH3_MONO; 
      }
      else{ // average
        mass += MASS_NH3_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, C_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for C ion");
        return FALSE;
      }
    }
    
    // add X ion
    if(constraint->ion_type == X_ION || constraint->ion_type == ALL_ION){
      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];

      if(constraint->mass_type == MONO){
        mass += MASS_CO_MONO + MASS_H2O_MONO;     
      }
      else{ // average
        mass += MASS_CO_AVERAGE + MASS_H2O_AVERAGE; 
      }
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, X_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge for X ion");
        return FALSE;
      }
    }
    
    // add Y ion
    if(constraint->ion_type == ALL_ION || 
       constraint->ion_type == BY_ION ||
       constraint->ion_type == BYA_ION ||
       constraint->ion_type == Y_ION){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint->mass_type == MONO){
        mass += MASS_H2O_MONO; 
      }
      else{ // average
        mass += MASS_H2O_AVERAGE;
      }
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, Y_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Y ion");
        return FALSE;
      }
      
    }
    
    // add Z ion
    if(constraint->ion_type == Z_ION ||
       constraint->ion_type == ALL_ION ){

      // set mass 
      mass = mass_matrix[(int)mass_matrix[0]] - mass_matrix[(int)mass_matrix[0] - cleavage_idx];
      
      if(constraint->mass_type == MONO){
        mass = mass - MASS_NH3_MONO + MASS_H2O_MONO; 
      }
      else{ // average
        mass = mass - MASS_NH3_AVERAGE + MASS_H2O_AVERAGE;
      }
      
      // add ions up to max charge
      if(!add_ions_by_charge(ion_series, mass, cleavage_idx, Z_ION)){
        carp(CARP_ERROR, "failed to add ions by different charge Z ion");
        return FALSE;
      }
    }
  }

  // add P ion(precursor ion)?
  if(ion_series->constraint->precursor_ion){
    
    // set mass 
    mass = mass_matrix[(int)mass_matrix[0]];
    
    // mass type
    if(constraint->mass_type == MONO){
      mass += MASS_H2O_MONO; 
    }
    else{ // average
      mass += MASS_H2O_AVERAGE;
    }
    
    // add ions up to max charge
    if(!add_ions_by_charge(ion_series, mass, (int)mass_matrix[0], P_ION)){
      carp(CARP_ERROR, "failed to add ions by different charge P ion");
      return FALSE;
    }
  }

  return TRUE;
}

/**
 * The modification depends on the loss/add && if the ion contains RKQN or STED
 * The number of losses possible cannot exceed the number of RKQN or STED in the ion
 * The loss_limit array in the ion_series must be populated prior to this method call
 *\returns TRUE if the ion can lose the mod_type modification, else FALSE
 */
BOOLEAN_T can_ion_lose_modification(
  ION_SERIES_T* ion_series, ///< ion_series to print -in/out  
  ION_T* ion, ///< the ion to check if can lose nh3 -in
  ION_MODIFICATION_T mod_type, ///< generate ions of this modification_type -in/out
  int increment  ///< the add/loss of the modification
  )
{
  int cleavage_idx = get_ion_cleavage_idx(ion);

  // check for NH3 modification
  if(mod_type == NH3){
    // adding is ok
    if(increment >= 0){
      return TRUE;
    }
    
    // is forward ion_type(ABC)?
    if(is_forward_ion_type(ion)){
      // does this ion contain enough counts of the RKQN
      if(-increment >  (&ion_series->loss_limit[cleavage_idx-1])->nh3){
        return FALSE;
      }
      return TRUE;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the RKQN
      if(cleavage_idx == ion_series->peptide_length){
        if(-increment > (&ion_series->loss_limit[ion_series->peptide_length-1])->nh3){
          return FALSE;
        }
      }
      else if(-increment >  
              ((&ion_series->loss_limit[ion_series->peptide_length-1])->nh3 - 
               (&ion_series->loss_limit[ion_series->peptide_length - cleavage_idx - 1])->nh3)){
        return FALSE;
      }
      return TRUE;
    }
  }

  // check for H2O modification
  if(mod_type == H2O){
    // adding is ok
    if(increment >= 0){
      return TRUE;
    }
    
    // is forward ion_type(ABC)?
    if(is_forward_ion_type(ion)){
      // does this ion contain enough counts of the STED
      if(-increment >  (&ion_series->loss_limit[cleavage_idx-1])->h2o){
        return FALSE;
      }
      return TRUE;
    }
    else{// backward ions XYZ
      // does this ion contain enough counts of the STED
      if(cleavage_idx == ion_series->peptide_length){
        if(-increment > (&ion_series->loss_limit[ion_series->peptide_length-1])->h2o){
          return FALSE;
        }
      }
      else if(-increment >  
              ((&ion_series->loss_limit[ion_series->peptide_length-1])->h2o - 
               (&ion_series->loss_limit[ion_series->peptide_length - cleavage_idx - 1])->h2o)){
        return FALSE;
      }
      return TRUE;
    }
  }
  
  // check for ISOTOPE modification
  else if(mod_type == ISOTOPE){
    // adding is ok
    if(increment >= 0 && get_ion_type(ion) != P_ION){
      return TRUE;
    }
    // add more constraint if needed
  }
  // check for FLANK modification
  else if(mod_type == FLANK){
    // only add flanking ions to type B,Y ions
    if(get_ion_type(ion) == B_ION || get_ion_type(ion) == Y_ION){
      // only add ions with no modifications
      if(!ion_is_modified(ion)){
        return TRUE;
      }
    }
  }
  
  return FALSE;
}

/**
 * creates all the ions with specific modifications up to the max charge
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T generate_ions(
  ION_SERIES_T* ion_series, ///< ion_series to print -in/out  
  ION_MODIFICATION_T mod_type ///< generate ions of this modification_type -in/out
  )
{
  int ion_idx = 0;
  int total_ion = ion_series->num_ions;
  ION_T* working_ion = NULL;
  ION_T* new_ion = NULL;
  int* modifications = ion_series->constraint->modifications;

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // check if there's enough space to add the new ions
  if(ion_series->num_ions*2 > MAX_IONS){
    carp(CARP_ERROR, "exceeds ion array size in ion_series");
    return FALSE;
  }

  // reset modification increment, if mod_type loss
  if(modifications[mod_type] < 0){
    type_increment = -1;
  }
  
  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ion_series->ions[ion_idx];
    
    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!(can_ion_lose_modification(ion_series, working_ion, mod_type, type_increment))){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; abs(type_idx) <= abs(modifications[mod_type]); ){
      // copy the src ion, into new ion
      new_ion = allocate_ion();
      copy_ion(working_ion, new_ion, get_ion_peptide_sequence(working_ion));
      
      // add the modification to the new ion
      add_modification(new_ion, mod_type, type_idx, ion_series->constraint->mass_type);
      
      // add ion to ion_series
      add_ion_to_ion_series(ion_series, new_ion);
     
      // can this ion generate a mod_type modification for the next count of modification?, 
      if(!(can_ion_lose_modification(ion_series, working_ion, mod_type, 
              (type_idx += type_increment)))){
        break;
      }
    }
  }
  return TRUE;
}

/**
 * creates all the flanking ions up to the max charge
 * can only create flanking ions that are B|Y ions and don't have modification
 * assumes the ions with no modification all are at the begining of the ion[] in ion_series
 * copies all the existing ions that can be modified,
 * then applies the different modifications then adds the new modified ions to ion_series
 *\returns TRUE if successfully generates all the ions with modifications, else FALSE
 */
BOOLEAN_T generate_ions_flank(
  ION_SERIES_T* ion_series ///< ion_series to print -in/out
  )
{
  int ion_idx = 0;
  int total_ion = ion_series->num_ions;
  ION_T* working_ion = NULL;
  ION_T* new_ion = NULL;
  ION_T* new_ion2 = NULL;
  int* modifications = ion_series->constraint->modifications;

  // modification index
  int type_idx = 0;
  int type_increment = 1;

  // check if there's enough space to add the new ions
  if(ion_series->num_ions*2 > MAX_IONS){
    carp(CARP_ERROR, "exceeds ion array size in ion_series");
    return FALSE;
  }
  
  // iterate over all the ions to determine which ones should be copies and modified
  for(; ion_idx < total_ion; ++ion_idx){
    working_ion = ion_series->ions[ion_idx];
    
    // no more ions that are not modified, thus done
    if(get_ion_single_modification_count(working_ion, NH3) != 0){
      break;
    }

    // can this ion generate a mod_type modification?, if not skip to next ion
    if(!can_ion_lose_modification(ion_series, working_ion, FLANK, type_increment)){      
      continue;
    }
     
    // add/sub thorugh all mod_type modifications!!!
    for(type_idx = type_increment; type_idx <= modifications[FLANK]; type_idx += type_increment){
      // copy the src ion, into new ion
      new_ion = allocate_ion();
      new_ion2 = allocate_ion();
      copy_ion(working_ion, new_ion, get_ion_peptide_sequence(working_ion));
      copy_ion(working_ion, new_ion2, get_ion_peptide_sequence(working_ion));
      
      // add the modification to the new ion
      add_modification(new_ion, FLANK, type_idx, ion_series->constraint->mass_type);
      add_modification(new_ion2, FLANK, -type_idx, ion_series->constraint->mass_type);

      // add ion to ion_series
      add_ion_to_ion_series(ion_series, new_ion);
      add_ion_to_ion_series(ion_series, new_ion2);
    }
  }
  return TRUE;
}

/**
 * \brief The engine of ion series. Predicts all the ions from the
 * peptide that meet the ion constraint. All predicted ions are stored
 * in the ion_series as ion objects. 
 */
void predict_ions(
  ION_SERIES_T* ion_series ///< the ion series to predict ions for -in
  )
{
  if(ion_series->is_predicted){
    carp(CARP_WARNING, "The ion series has already been predicted and added");
    return;
  }

  ION_CONSTRAINT_T* constraint = ion_series->constraint;
  
  // create a mass matrix
  FLOAT_T* mass_matrix = 
    create_ion_mass_matrix(ion_series->modified_aa_seq, constraint->mass_type, ion_series->peptide_length);  
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
  scan_for_aa_for_neutral_loss(ion_series);
  
  // generate ions without any modifications
  if(!generate_ions_no_modification(ion_series, mass_matrix)){
    carp(CARP_FATAL, "failed to generate ions, no modifications");
  }

  // create modification ions?
  if(ion_series->constraint->use_neutral_losses){
    
    // generate ions with nh3 modification
    if(abs(constraint->modifications[NH3]) > 0){
      if(!generate_ions(ion_series, NH3)){
        carp(CARP_FATAL, "failed to generate ions, NH3 modifications");
      }
    }
    
    // generate ions with h2o modification
    if(abs(constraint->modifications[H2O]) > 0){
      if(!generate_ions(ion_series, H2O)){
        carp(CARP_FATAL, "failed to generate ions, H2O modifications");
      }
    }

    // generate ions with isotope modification
    if(constraint->modifications[ISOTOPE] > 0){
      if(!generate_ions(ion_series, ISOTOPE)){
        carp(CARP_FATAL, "failed to generate ions, ISOTOPE modifications");
      }
    }

    // generate ions with flank modification
    if(constraint->modifications[FLANK] > 0){
      if(!generate_ions_flank(ion_series)){
        carp(CARP_FATAL, "failed to generate ions, FLANK modifications");
      }
    }
    
    // add more modifications here

  }
  
  // ion series now been predicted
  ion_series->is_predicted = TRUE;

  // free mass matrix
  free(mass_matrix);
}
/**
 * Assign peaks to the nearest ions, within a tolerance (set in param file)
 */
void ion_series_assign_nearest_peaks(
    ION_SERIES_T* ion_series, 
    SPECTRUM_T* spectrum){

  //  FLOAT_T max = 0.5; // TODO set in param file 
  FLOAT_T max = get_double_parameter("ion-tolerance"); 
  ION_T* ion = NULL;
  ION_ITERATOR_T* iterator = new_ion_iterator(ion_series);
  PEAK_T* peak = NULL;
  while(ion_iterator_has_next(iterator)){
    ion = ion_iterator_next(iterator);
    FLOAT_T mz = get_ion_mass_z(ion); // TODO change to mz, not mass_z
    peak = get_nearest_peak(spectrum, mz, max);
    set_ion_peak(ion, peak);
  }
  free_ion_iterator(iterator);
}

/**
 * Copies ion_series object from src to dest.
 *  must pass in a memory allocated ION_SERIES_T* dest
 * does not copy the loss_limit.
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
  dest->peptide_length = src->peptide_length;
  //mod seq???

  // add copy of pointer ion constraint
  dest->constraint = src->constraint;

  // add copy ion, add ion_filtered_iterator
  ION_ITERATOR_T* iterator = new_ion_iterator(src);

  // iterate over all ions in src and copy them into dest
  while(ion_iterator_has_next(iterator)){
    src_ion = ion_iterator_next(iterator);
    // add ion
    dest_ion = allocate_ion(); 
    copy_ion(src_ion, dest_ion, dest->peptide);
    add_ion_to_ion_series(dest, dest_ion);
  }
  // free up iterator
  free_ion_iterator(iterator);

  dest->is_predicted = TRUE;
}

/*************************************
 * ION_SERIES_T: get and set methods
 ************************************/

/**
 * \returns the ion that meets the criteria or NULL
 * TODO possibly reimplement if linear scan is too slow
 */
ION_T* get_ion_series_ion(
    ION_SERIES_T* ion_series,
    ION_CONSTRAINT_T* ion_constraint,
    int cleavage_idx
    ){

  ION_FILTERED_ITERATOR_T* ion_iterator = 
    new_ion_filtered_iterator(ion_series, ion_constraint);
  ION_T* ion = NULL;

  while(ion_filtered_iterator_has_next(ion_iterator) == TRUE){
    ion = ion_filtered_iterator_next(ion_iterator);
    if(get_ion_cleavage_idx(ion) == cleavage_idx){
      return ion;
    }
  }
  return NULL;
}

/**
 *\returns the peptide length of which the ions are made
 */
int get_ion_series_peptide_length(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  )
{

  return ion_series->peptide_length;
}

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
  // free previous sequence
  if(ion_series->peptide != NULL){
    free(ion_series->peptide);
  }

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
 * set the of the ion_series object
 */
void set_ion_series_ion_constraint(
  ION_SERIES_T* ion_series, ///< the working ion_series -in
  ION_CONSTRAINT_T* constraint///<  -in
  )
{
  ion_series->constraint = constraint;
}

/**
 *\returns the total number of ions in the ion_series object
 */
int get_ion_series_num_ions(
  ION_SERIES_T* ion_series ///< the working ion_series -in                          
  )
{
  return ion_series->num_ions;
}

/**
 *\returns the total number of ion_type in the ion_series object
 */
int get_ion_series_num_ions_one_type(
  ION_SERIES_T* ion_series, ///< the working ion_series -in                          
  ION_TYPE_T ion_type ///< the type of ions -in
  )
{
  return ion_series->num_specific_ions[ion_type];
}

/*************************
 * ION_CONSTRAINT methods
 *************************/

/**
 *\returns an empty heap allocated ion_constraint
 */
ION_CONSTRAINT_T* allocate_ion_constraint(void){
  ION_CONSTRAINT_T* constraint = (ION_CONSTRAINT_T*)mycalloc(1, sizeof(ION_CONSTRAINT_T));
  int modification_idx = 0;

  // initialize all modifications count to 0
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    constraint->modifications[modification_idx] = 0;
  }
  return constraint;
}

/**
 * modification, all modifications 0
 * add more modifications as needed using the set_ion_constraint_modification
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< max charge of the ions <= the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion  ///< should include precursor ion?
  )
{
  ION_CONSTRAINT_T* constraint = allocate_ion_constraint();
  constraint->use_neutral_losses = FALSE;
  
  // set all fields of constraint
  constraint->mass_type = mass_type;
  constraint->max_charge = max_charge;
  constraint->min_charge = 0;
  constraint->exact_modifications = FALSE;
  constraint->ion_type = ion_type;
  constraint->precursor_ion = precursor_ion;
  constraint->pointer_count = 1;

  return constraint;
}

/**
 * \brief Create a new ion constraint based on the score type and the
 * charge of the peptide to be modeled.  Uses other
 * new_ion_constraint_ methods for some types.
 *
 * \returns A newly allocated ion constraint.
 */
ION_CONSTRAINT_T* new_ion_constraint_smart(
  SCORER_TYPE_T score_type,
  int charge
){
  ION_CONSTRAINT_T* new_constraint = NULL;

  switch(score_type){
  case SP:
    new_constraint = new_ion_constraint_sequest_sp(charge);
    break;
  case XCORR:
    new_constraint = new_ion_constraint_sequest_xcorr(charge);
    break;
  default:
    // use default type for others
    new_constraint = 
      new_ion_constraint(get_mass_type_parameter("fragment-mass"),
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
ION_CONSTRAINT_T* new_ion_constraint_gmtk(
  int charge 
  )
{
  ION_CONSTRAINT_T* constraint = NULL;

  int max_charge = 1;
  if(charge > 1){
    max_charge = charge;
  }  
  constraint = new_ion_constraint(MONO, max_charge, ALL_ION, FALSE);

  // set all modifications count for gmtk
  constraint->use_neutral_losses = TRUE;
  constraint->min_charge = 1;
  constraint->modifications[NH3] = -1;
  constraint->modifications[H2O] = -1;
  constraint->modifications[ISOTOPE] = 0;
  constraint->modifications[FLANK] = 0;

  return constraint;
}


/**
 * modification, sets all fields for sequest settings
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest(
  MASS_TYPE_T mass_type, ///< the mass_type to use MONO|AVERAGE
  int max_charge, ///< the maximum charge of the ions. 
                  ///< cannot exceed the parent peptide's charge
  ION_TYPE_T ion_type, ///< the ion types the peptide series should include
  BOOLEAN_T precursor_ion ///< should include precursor ion?
  )
{
  ION_CONSTRAINT_T* constraint = new_ion_constraint(mass_type, max_charge, ion_type, precursor_ion);

  // set                                                     
  constraint->use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  constraint->modifications[NH3] = 1;
  constraint->modifications[H2O] = 1;
  constraint->modifications[ISOTOPE] = 1;
  constraint->modifications[FLANK] = 1;
  
  return constraint;
}

/**
 * modification, sets all fields for sequest Sp scoring settings
 * make B, Y type ions
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest_sp(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  ION_CONSTRAINT_T* constraint = NULL;
  // charge = 1;
  if(charge == 1){
    constraint = new_ion_constraint(MONO, 1, BY_ION, FALSE);
  }  
  else{
    --charge;
    constraint = new_ion_constraint(MONO, charge, BY_ION, FALSE);
  }
  
  // set                                                     
  constraint->use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  constraint->modifications[NH3] = 0;
  constraint->modifications[H2O] = 0;
  constraint->modifications[ISOTOPE] = 0;
  constraint->modifications[FLANK] = 0;
  
  return constraint;
}


/**
 * modification, sets all fields for Sequest Xcorr scoring settings
 * make B, Y, A type ions
 *\returns a new heap allocated ion_constraint
 */
ION_CONSTRAINT_T* new_ion_constraint_sequest_xcorr(
  int charge ///< the maximum charge of the ions, cannot exceed the parent peptide's charge
  )
{
  ION_CONSTRAINT_T* constraint = NULL;
  // charge = 1;
  if(charge == 1){
    constraint = new_ion_constraint(MONO, 1, BYA_ION, FALSE);
  }  
  else{
    --charge;
    constraint = new_ion_constraint(MONO, charge, BYA_ION, FALSE);
  }
  
  // set                                                     
  constraint->use_neutral_losses = TRUE;
  
  // set all modifications count for sequest
  constraint->modifications[NH3] = 0;// -1;
  constraint->modifications[H2O] = 0;// -1;
  constraint->modifications[ISOTOPE] = 0;// not sure
  constraint->modifications[FLANK] = 0;
  
  return constraint;
}

/**
 * Frees an allocated ion_constraint object.
 */
void free_ion_constraint(
  ION_CONSTRAINT_T* ion_constraint///< the ion constraints to enforce -in
  )
{
  ion_constraint->pointer_count--;
  if (ion_constraint->pointer_count == 0){
    free(ion_constraint);
  }
}


/**
 * copies ion_constraint pointer
 */
ION_CONSTRAINT_T* copy_ion_constraint_ptr(
    ION_CONSTRAINT_T* ion_constraint
    ){
  ion_constraint->pointer_count++;
  return ion_constraint;
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

  // if use natural loss, copy
  if(src->use_neutral_losses){
    // iterate over all modifications a update new constraint
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      dest->modifications[modification_idx] = src->modifications[modification_idx];
    }
  }
  
  dest->mass_type = src->mass_type;
  dest->max_charge = src->max_charge;
  dest->ion_type = src->ion_type;
  dest->precursor_ion = src->precursor_ion;
}

/** FIXME!!!! double check
 * Determines if a ion satisfies a ion_constraint.
 * \returns TRUE if the constraint is satisified. FALSE if not.
 */
BOOLEAN_T ion_constraint_is_satisfied(
   ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
   ION_T* ion ///< query ion -in
   )
{
  int* counts = NULL;

  // TODO Fix
  BOOLEAN_T return_val = TRUE;
  // print_ion(ion, stderr);
  // fprintf(stderr, "%i->%i\n", ion_constraint->min_charge, ion_constraint->max_charge);
  // check ion type
  ION_TYPE_T ion_type = get_ion_type(ion);
  if(
     !(ion_type == ion_constraint->ion_type)
      
     &&
     
     !((ion_constraint->ion_type == BY_ION) && (ion_type == B_ION || ion_type == Y_ION)) 
     
     &&

     !((ion_constraint->ion_type == BYA_ION) 
          && 
       (ion_type == B_ION || ion_type == Y_ION || ion_type == A_ION)) 
     
     &&
     
     !(ion_constraint->ion_type == ALL_ION)

     ){
     
    // precursor ion?
    if(!(ion_constraint->precursor_ion && ion_type == P_ION)){
      return_val = FALSE;
    }
  }
  
  // check charge
  if(get_ion_charge(ion) > ion_constraint->max_charge){
    return_val = FALSE;
  }
  
  if(get_ion_charge(ion) < ion_constraint->min_charge){
    return_val = FALSE;
  }

  // check modifications
  counts = get_ion_modification_counts(ion);
  int mod_idx;
  for(mod_idx=0; mod_idx < MAX_MODIFICATIONS; ++mod_idx){
    if(ion_constraint->modifications[mod_idx] >= 0){
      if(counts[mod_idx] > ion_constraint->modifications[mod_idx]){
        return_val = FALSE;
        break;
      }
    }
    else{
      if(counts[mod_idx] < ion_constraint->modifications[mod_idx]){
        return_val = FALSE;
        break;
      }
    }
    if (ion_constraint->exact_modifications){
      if(counts[mod_idx] != ion_constraint->modifications[mod_idx]){
        return_val = FALSE; 
        break;
      }
    }
  }
  
  // Add more checks here as more contraints are added

  // fprintf(stderr, "r = %i\n", return_val);
  return return_val;
}


/**
 * sets the modification count
 * can only add isotopes
 */
void set_ion_constraint_modification(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  ION_MODIFICATION_T mod_type, ///< ion modification type -in
  int count  ///< the count of the modification -in  
  )
{
  // set modification count, can only add for isotope
  if(mod_type != ISOTOPE){
    ion_constraint->modifications[mod_type] = count;
  }
  else{
    ion_constraint->modifications[mod_type] = abs(count);
  }

  // set neutral loss to TRUE if needed
  if(!ion_constraint->use_neutral_losses){
    ion_constraint->use_neutral_losses = TRUE;
  }
}

/**
 * sets the exact modification boolean 
 */
void set_ion_constraint_exactness(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  BOOLEAN_T exactness ///< whether to use exact mods or not -in
  ){
  ion_constraint->exact_modifications = exactness;
  if (exactness == TRUE){
    ion_constraint->min_charge = ion_constraint->max_charge;
  }
}
 
/**
 * gets the modification count for specific mod_type
 */
int get_ion_constraint_modification(
  ION_CONSTRAINT_T* ion_constraint,///< the ion constraints to enforce -in
  ION_MODIFICATION_T mod_type ///< ion modification type -in
  )
{
  return ion_constraint->modifications[mod_type];
}

/**
 * gets the mass type of the ion_constraint
 */
MASS_TYPE_T get_ion_constraint_mass_type(
  ION_CONSTRAINT_T* ion_constraint ///< the ion constraints to enforce -in
  )
{
  return ion_constraint->mass_type;
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
  return ion_iterator->ion_series->num_ions > ion_iterator->ion_idx;
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

  // iterate over ions until discovers the first ion that meets the ion constraint
  while(ion_iterator->ion_idx < ion_iterator->array_size){
    // get next ion
    ion = ion_iterator->ion_array[ion_iterator->ion_idx];
    
    // check if the current ion satisfies the ion_constraint for the iterator
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
  
  // set constraint, ion_series
  iterator->constraint = constraint;
  iterator->ion_series = ion_series;
  iterator->has_next = FALSE;

  // set the working array of ions
  if(constraint->ion_type == ALL_ION ||
     constraint->ion_type == BY_ION ||
     constraint->ion_type == BYA_ION){

    iterator->ion_array = ion_series->ions;
    iterator->array_size = ion_series->num_ions;
  }
  else{
    iterator->ion_array = ion_series->specific_ions[constraint->ion_type];
    iterator->array_size = ion_series->num_specific_ions[constraint->ion_type];
  }
  
  // initialize iterator
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
  
  // check if a ion is present to return
  if(!ion_iterator->has_next){
    carp(CARP_FATAL, "index out of bounds for ion_filtered_iterator");
  }
  
  next_ion = ion_iterator->ion;
  
  // re-initialize iterator
  setup_ion_filtered_iterator(ion_iterator);

  return next_ion;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
