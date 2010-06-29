/*************************************************************************//**
 * \file ion.cpp
 * \brief Object for representing a single ion.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "objects.h"
#include "ion.h"
#include "alphabet.h"
#include "peptide.h"
#include "peak.h"
#include "mass.h"
#include "utils.h"
#include <sys/types.h>
#include <netinet/in.h>
#include <inttypes.h>

static const int MZ_INT_MAX = 10;
static const int MZ_INT_MIN = 0;
static const int PAIRED_ION_INTS = 12;
static const int PAIRED_ION_FLOATS = 6;
static const int SINGLE_ION_FLOATS = 3;
static const int SINGLE_ION_INTS = 9;

// At one point I need to reverse the endianness for pfile_create to work
// Apparently that is no longer true. Hence 0 below.
static const bool REVERSE_ENDIAN = 0;

/**
 * Array to store the modification masses
 */
FLOAT_T modification_masses[MAX_MODIFICATIONS];

/**
 * Have we initialized the modification_masses?
 */
BOOLEAN_T initialized_modification_masses = FALSE;

static const int DETECTABLE_MZ_MIN = 200;
static const int DETECTABLE_MZ_MAX = 2400;

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
  FLOAT_T peptide_mass; ///< the mass of the peptide. For efficiency
  int modification_counts[MAX_MODIFICATIONS]; ///< an array of the number of different ion modifications
  FLOAT_T ion_mass_z;   ///< The mass/z of the ion. 
  PEAK_T* peak;  ///< The assigned peak. NULL if no peak // TODO add ptr count?
};

/**
 * initializes the mass array
 */
void initialize_modification_masses(
   MASS_TYPE_T mass_type ///< mass type (average, mono) -in
)
{
  // set modification mass
  if(mass_type == MONO){
    modification_masses[NH3] = MASS_NH3_MONO;
    modification_masses[H2O] = MASS_H2O_MONO ;
    modification_masses[ISOTOPE] = 1; // FIXME check this!!!
    modification_masses[FLANK] = 1; // FIXME check this!!!
  }
  else if(mass_type == AVERAGE){
    modification_masses[NH3] = MASS_NH3_AVERAGE;
    modification_masses[H2O] = MASS_H2O_AVERAGE;
    modification_masses[ISOTOPE] = 1; // FIXME check this!!!
    modification_masses[FLANK] = 1; // FIXME check this!!!
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
  // allocate new ion
  ION_T* ion = (ION_T*)mycalloc(1, sizeof(ION_T));
  ion->type = type;
  ion->cleavage_idx = cleavage_idx;
  ion->charge = charge;
  ion->peptide_sequence = peptide;
  // TODO get mass type from param file
  ion->peptide_mass = calc_sequence_mass(peptide, MONO); 
  ion->peak = NULL;
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
  // get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  // calculate and set ion mass/z
  if(!calc_ion_mass_z(ion, mass_type, FALSE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
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
  // get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = modification_counts[modification_idx];
  }
  
  // calculate and set ion mass/z
  if(!calc_ion_mass_z(ion, mass_type, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
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
  FLOAT_T base_mass, ///< the base mass of the ion -in
  int* modification_counts ///< an array of modification counts for each modification -in
  )
{
  // get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = modification_counts[modification_idx];
  }
  
  // calculate and set ion mass/z
  if(!calc_ion_mass_z_with_mass(ion, mass_type, base_mass, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
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
  FLOAT_T base_mass ///< the base mass of the ion -in
  )
{
  // get new basic ion
  ION_T* ion = new_basic_ion(type, cleavage_idx, charge, peptide);
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    ion->modification_counts[modification_idx] = 0;
  }
  
  // calculate and set ion mass/z
  if(!calc_ion_mass_z_with_mass(ion, mass_type, base_mass, TRUE)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
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
  // print all fields of ion
  fprintf(file, "%.2f\t%.2f\t%d\t%d\t%d", 
      ion->ion_mass_z, (ion->ion_mass_z)*ion->charge, ion->charge, 
      (int)ion->type, ion->cleavage_idx);

  // iterate over all modification counts
  int mod_idx;
  for(mod_idx=0; mod_idx < MAX_MODIFICATIONS; ++mod_idx){
    fprintf(file,"\t%d", ion->modification_counts[mod_idx]);
  }
  
  fprintf(file,"\n");
}


/**
 * prints the ION_T object to the file, in the
 * following format for GMTK single-ion models:
 *
 * 1. m/z ratio - ratio of the ion's mass-to-charge to the peptide's m/z
 * 2. raw - raw intensity
 * 3. rank - the ion rank
 * 4. proton mobility - always set to 1 (for now) FIXME
 * 5. m/z ratio int
 * 6. index of the amide bond cleavage from N-term
 * 7. index of the amide bond cleavage from C-term
 * 8. Left amino acid ID
 * 9. Right amino acid ID
 * 10. Is this ion possible?
 * 11. Is this ion detectable?
 * 12. Is this ion detected?
 */
void print_ion_gmtk_single(
  ION_T* ion, ///< print this ion -in
  FILE* file  ///< to this file -in
  ){

  int has_mobile_proton = 1;
  int is_detectable = 0;
  int is_detected = 0;
  int is_possible = 0;

  FLOAT_T intensity = 0.0;
  FLOAT_T intensity_rank = 0.0;
  if (ion->peak != NULL){
    intensity = get_peak_intensity(ion->peak);
    intensity_rank = get_peak_intensity_rank(ion->peak);
    is_detected = 1;
  }

  if ((ion->ion_mass_z >= DETECTABLE_MZ_MIN) &&
      (ion->ion_mass_z <= DETECTABLE_MZ_MAX)){
    is_detectable = 1;
  }

  FLOAT_T mz_ratio = (ion->ion_mass_z)/(ion->peptide_mass);
  int mz_int = (int)(mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);

  const char* format = "%.6f\t%.6f\t%.6f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n";
  fprintf(file, format,
      mz_ratio,                                                 // 1 
      intensity,                                                // 2 
      intensity_rank,                                           // 3 
      has_mobile_proton,                                        // 4 
      mz_int,                                                   // 5 
      ion->cleavage_idx,                                        // 6
      strlen(ion->peptide_sequence) - ion->cleavage_idx + 1,    // 7
      amino_to_int(ion->peptide_sequence[ion->cleavage_idx-1]), // 8 
      amino_to_int(ion->peptide_sequence[ion->cleavage_idx]),   // 9 
      is_possible,                                              // 10 
      is_detectable,                                            // 11
      is_detected                                               // 12
  );
}

/*
 * Same output as above, but in binary
 *
 * ints 
 *
 * 0. m/z ratio - ratio of the ion's mass-to-charge to the peptide's m/z
 * 1. raw - raw intensity
 * 2. rank - the ion rank
 *
 * floats 
 *
 * 0. proton mobility - always set to 1 (for now) FIXME
 * 1. m/z ratio int
 * 2. index of the amide bond cleavage from N-term
 * 3. index of the amide bond cleavage from C-term
 * 4. Left amino acid ID 
 * 5. Right amino acid ID
 * 6. Is this ion possible?
 * 7. Is this ion detectable?
 * 8. Is this ion detected?

 */
void print_ion_gmtk_single_binary(
  ION_T* ion, ///< print this ion -in
  FILE* file,  ///< to this file -in
  int sentence_idx,
  int frame_idx
  ){

  FLOAT_T* float_array = (FLOAT_T*)mycalloc(sizeof(FLOAT_T), SINGLE_ION_FLOATS);
  int* int_array = (int*)mycalloc(sizeof(int), SINGLE_ION_INTS);

  FLOAT_T mz_ratio = (ion->ion_mass_z)/(ion->peptide_mass);
  float_array[0] = mz_ratio;                              // 0
  float_array[1] = 0.0;                                   // 1
  float_array[2] = 0.0;                                   // 2
  int is_detected = 0;
  if (ion->peak != NULL){
    float_array[1] = get_peak_intensity(ion->peak);       // 1 
    float_array[2] = get_peak_intensity_rank(ion->peak);  // 2 
    is_detected = 1; 

#ifdef LOG_INTENSITY
    // START
    // replace the rank intensity with the negative log of the 
    // TIC normalized intensity. A hack, I know
    float_array[2] = 10.0 + log(float_array[1]);
#endif

  }

  int mz_int = (int)(mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int cterm_idx = strlen(ion->peptide_sequence) - ion->cleavage_idx; 
  int left_amino = amino_to_int(ion->peptide_sequence[ion->cleavage_idx-1]);
  int right_amino = amino_to_int(ion->peptide_sequence[ion->cleavage_idx]);
  int is_detectable = 0;
  if ( 
       ((ion->ion_mass_z >= DETECTABLE_MZ_MIN) 
              &&
        (ion->ion_mass_z <= DETECTABLE_MZ_MAX)) 

         || 
         
       is_detected

     ){
    is_detectable = 1;                                    
  }


  int_array[0] = 1;                 // 0
  int_array[1] = mz_int;            // 1
  int_array[2] = ion->cleavage_idx; // 2
  int_array[3] = cterm_idx;         // 3
  int_array[4] = left_amino;        // 4
  int_array[5] = right_amino;       // 5
  int_array[6] = 1;                 // 6
  int_array[7] = is_detectable;     // 7
  int_array[8] = is_detected;       // 8

  if (REVERSE_ENDIAN){

    int idx;
    for (idx=0; idx < 9; idx++){
      int_array[idx] = htonl(int_array[idx]);
    }
    for (idx=0; idx < 3; idx++){
      // htonl does not seem to work on the floats (!!)
      // so I will reverse the bytes by hand
      FLOAT_T old_float;
      FLOAT_T new_float;

      old_float = float_array[idx];

      char *forward_endian = (char*) &old_float;
      char *reversed_endian = (char*) &new_float;

      reversed_endian[0] = forward_endian[3];
      reversed_endian[1] = forward_endian[2];
      reversed_endian[2] = forward_endian[1];
      reversed_endian[3] = forward_endian[0];
      float_array[idx] = new_float;
    }
    
    int big_end_sentence_idx = htonl(sentence_idx);
    int big_end_frame_idx = htonl(frame_idx);

    fwrite(&big_end_sentence_idx, sizeof(int), 1, file);
    fwrite(&big_end_frame_idx, sizeof(int), 1, file);

  } else {

    fwrite(&sentence_idx, sizeof(int), 1, file);
    fwrite(&frame_idx, sizeof(int), 1, file);

  }

  fwrite(float_array, sizeof(FLOAT_T), 3, file);
  fwrite(int_array, sizeof(int), 9, file);
  free(float_array);
  free(int_array);
}

/**
 * A hack routine to print out a null ion if there are none in the series.
 * For using neutral losses with GMTK.
 * Come in both binary and ascii versions.
 */
void print_null_ion_gmtk_single_binary(
  FILE* file,
  int sentence_idx,
  int frame_idx
  ){

  FLOAT_T float_array[3] = {0.0, 0.0, 0.0};
  int int_array[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  fwrite(&sentence_idx, sizeof(int), 1, file);
  fwrite(&frame_idx, sizeof(int), 1, file);
  fwrite(float_array, sizeof(FLOAT_T), 3, file);
  fwrite(int_array, sizeof(int), 9, file);
}

/**
 * A hack routine to print out a null ion if there are none in the series.
 * For using neutral losses with GMTK.
 */
void print_null_ion_gmtk_paired_binary(
  FILE* file,
  int sentence_idx,
  int frame_idx
  ){

  // FIX get rid of magic numbers
  FLOAT_T float_array[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  int int_array[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  fwrite(&sentence_idx, sizeof(int), 1, file);
  fwrite(&frame_idx, sizeof(int), 1, file);
  fwrite(float_array, sizeof(FLOAT_T), 6, file);
  fwrite(int_array, sizeof(int), 12, file);
}


/**
 * A hack routine to print out a null ion if there are none in the series.
 * For using neutral losses with GMTK.
 * Come in both binary and ascii versions.
 */
void print_null_ion_gmtk_single(
  FILE* file
  ){

  const char* string = "0.0\t0.0\t0.0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
  fprintf(file, "%s", string);
}


/**
 * prints the location and fields of the two ION_T objects to the file
 * following format for GMTK paired-ion models:
 *
 * ints 
 *
 * 1.  m/z ratio int (from N-term)
 * 2.  m/z ratio int (from C-term)
 * 3.  peptide idx (from N-term)
 * 4.  peptide idx (from C-term)
 * 5.  aa Id (N-term)
 * 6.  aa Id (C-term)
 * 7.  first possible
 * 8.  first detectable
 * 9.  first detected
 * 10. second possible
 * 11. second observable
 * 12. second detected
 *
 * floats
 *
 * 1. m/z ratio FLOAT_T (from N-term)
 * 2. m/z ratio FLOAT_T (from C-term)
 * 3. first raw
 * 4. second raw
 * 5. first rank
 * 6. second rank
 *
 */
void print_ion_gmtk_paired_binary(
  ION_T* first_ion,   ///< print this ion -in
  ION_T* second_ion,  ///< print this ion -in
  FILE* file,         ///< to this file -in
  int sentence_idx, 
  int frame_idx
  ){
  
  FLOAT_T* float_array = (FLOAT_T*)mycalloc(sizeof(FLOAT_T), PAIRED_ION_FLOATS);
  int* int_array = (int*)mycalloc(sizeof(int), PAIRED_ION_INTS);

  // start with the floats
  FLOAT_T n_mz_ratio = (first_ion->ion_mass_z)/(first_ion->peptide_mass);
  float_array[0] = n_mz_ratio;                                    // 0
  // TODO 
  // subtract from 1.0?
  FLOAT_T c_mz_ratio =  1.0 - n_mz_ratio;                           // 1
  float_array[1] = c_mz_ratio;

  int first_is_detected = 0;
  if (first_ion->peak != NULL){
    // put in LOG_INTENSITY ?
    float_array[2] = get_peak_intensity(first_ion->peak);         // 2 
    float_array[4] = get_peak_intensity_rank(first_ion->peak);    // 4 
    first_is_detected = 1; 
  }

  int second_is_detected = 0;
  if (second_ion->peak != NULL){
    // put in LOG_INTENSITY ?
    float_array[3] = get_peak_intensity(second_ion->peak);        // 3 
    float_array[5] = get_peak_intensity_rank(second_ion->peak);   // 5 
    second_is_detected = 1; 
  }

  // next do the ints
  int n_mz_int = (int)(n_mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int c_mz_int = (int)(c_mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int cterm_idx = strlen(first_ion->peptide_sequence) 
    - first_ion->cleavage_idx; 
  int left_amino = amino_to_int(
      first_ion->peptide_sequence[first_ion->cleavage_idx-1]);
  int right_amino = amino_to_int(
      first_ion->peptide_sequence[first_ion->cleavage_idx]);
  int first_is_detectable = 0;
  if ( 
       ((first_ion->ion_mass_z >= DETECTABLE_MZ_MIN) 
              &&
        (first_ion->ion_mass_z <= DETECTABLE_MZ_MAX)) 

         || 
         
       first_is_detected

     ){
    first_is_detectable = 1;                                    
  }

  int second_is_detectable = 0;
  if ( 
       ((second_ion->ion_mass_z >= DETECTABLE_MZ_MIN) 
              &&
        (second_ion->ion_mass_z <= DETECTABLE_MZ_MAX)) 

         || 
         
       second_is_detected

     ){
    second_is_detectable = 1;                                    
  }

  int_array[0] = n_mz_int; 
  int_array[1] = c_mz_int; 
  int_array[2] = first_ion->cleavage_idx;
  int_array[3] = cterm_idx;
  int_array[4] = left_amino;
  int_array[5] = right_amino;
  int_array[6] = 1;  // TODO is possible always true?
  int_array[7] = first_is_detectable; 
  int_array[8] = first_is_detected; 
  int_array[9] = 1;  // TODO is possible always true?
  int_array[10] = second_is_detectable; 
  int_array[11] = second_is_detected; 

  // account for endian differences
  if (REVERSE_ENDIAN){

    int idx;
    for (idx=0; idx < PAIRED_ION_INTS; idx++){
      int_array[idx] = htonl(int_array[idx]);
    }
    for (idx=0; idx < PAIRED_ION_FLOATS; idx++){
      // htonl does not seem to work on the floats (!!)
      // so I will reverse the bytes by hand
      FLOAT_T old_float;
      FLOAT_T new_float;

      old_float = float_array[idx];

      char *forward_endian = (char*) &old_float;
      char *reversed_endian = (char*) &new_float;

      reversed_endian[0] = forward_endian[3];
      reversed_endian[1] = forward_endian[2];
      reversed_endian[2] = forward_endian[1];
      reversed_endian[3] = forward_endian[0];
      float_array[idx] = new_float;
    }
    
    sentence_idx = htonl(sentence_idx);
    frame_idx = htonl(frame_idx);
  } 

  // output to the file
  fwrite(&sentence_idx, sizeof(int), 1, file);
  fwrite(&frame_idx, sizeof(int), 1, file);
  fwrite(float_array, sizeof(FLOAT_T), PAIRED_ION_FLOATS, file);
  fwrite(int_array, sizeof(int), PAIRED_ION_INTS, file);
  free(float_array);
  free(int_array);
};

/**
 *\return the modified mass_z accodring to the modification type
 */
FLOAT_T modify_ion_mass_z(
  FLOAT_T mass_z, ///< the pre-modified ion mass_z -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  int charge, ///< charge of the ion
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // initalize the modified masses(average|mono);
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
  // update modification count
  ion->modification_counts[(int)modification] += modification_count;  
  // reset ion mass_z
  ion->ion_mass_z =  modify_ion_mass_z(ion->ion_mass_z, 
                                       modification_count, 
                                       modification, 
                                       ion->charge, 
                                       mass_type);
}

/**
 * Adds the given ION_MODIFICATION to this ion
 */
void set_ion_peak(
  ION_T* ion, ///< ion to which to add the peak -mod
  PEAK_T* peak ///< peak to add to this ion -in
  ){
  ion->peak = peak;
};

/**
 *\returns the ion's AA mass added all up
 */
FLOAT_T get_ion_mass(
  ION_T* ion, ///< the working ion -out/in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )                   
{
  char* ion_sequence = NULL;
  ION_TYPE_T ion_type = ion->type;
  int ion_length = 0;
  BOOLEAN_T memory_used = FALSE;
  FLOAT_T mass = 0;
  BOOLEAN_T reverse = FALSE;

  // get sequence for x,y,z ion
  if(ion_type == X_ION ||ion_type == Y_ION || ion_type == Z_ION){
    // convert the cleavage index into the actually index that start from the left.
    int real_cleavage_idx = strlen(ion->peptide_sequence) - ion->cleavage_idx;
    ion_sequence = &(ion->peptide_sequence[real_cleavage_idx]);
    ion_length = strlen(ion_sequence);
    reverse = TRUE;
  }
  // get sequence for a,b,c ion
  else{
    ion_length = ion->cleavage_idx;
    ion_sequence = (char*)mycalloc(ion_length+1, sizeof(char));
    strncpy(ion_sequence, ion->peptide_sequence, ion_length);
    memory_used = TRUE;
  }

  // add up all AA mass
  int ion_idx = 0;
  for(; ion_idx < ion_length; ++ion_idx){
    mass += get_mass_amino_acid(ion_sequence[ion_idx], mass_type);
  }
  
  // free ion sequence, only if memory allocated (a,b,c ions);
  if(memory_used){
    free(ion_sequence);
  }

  // if X,Y,Z ion add H2O
  if(reverse){
    if(mass_type == AVERAGE){
      mass += MASS_H2O_AVERAGE;
      //return mass + MASS_H2O_AVERAGE;
    }else{
      mass += MASS_H2O_MONO;
      //return mass + MASS_H2O_MONO;
    }
  }
  
  return mass;
}

/**
 *\return the modified mass accodring to the modification type
 */
FLOAT_T modify_ion_mass(
  FLOAT_T mass, ///< the pre-modified ion mass -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // initalize the modified masses(average|mono);
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
  FLOAT_T mass, ///< the basic mass of the ion -in
  BOOLEAN_T is_modified ///< are there any modifications for this ion? -in
  )
{
  FLOAT_T h_mass = MASS_H_MONO;

  // alter mass according to the modification
  if(is_modified){
    // iterate over all type of modifications for ion
    int modification_idx = 0;
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      // update ion mass if modification is needed
      if(ion->modification_counts[modification_idx] != 0){
        mass = modify_ion_mass(mass, ion->modification_counts[modification_idx], 
                               (ION_MODIFICATION_T)modification_idx, mass_type);
      }
    }
  }
  
  // reset h mass if needed
  if(mass_type == AVERAGE){
    h_mass = MASS_H_AVERAGE;
  }

  // convert mass to m/z and assigned to ion
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
  FLOAT_T mass = get_ion_mass(ion, mass_type);
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
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    dest->modification_counts[modification_idx] = src->modification_counts[modification_idx];
  }
  
  dest->ion_mass_z = src->ion_mass_z;
  dest->peptide_sequence = peptide_sequence;
  dest->peptide_mass = src->peptide_mass;
}


/**
 * \returns TRUE if forward ion_type(A,B,C), 
 * else reverse ion_type(X,Y,Z) FALSE
 */
BOOLEAN_T is_forward_ion_type(
  ION_T* ion ///< the ion to check if can lose nh3 -in                         
  )
{
  // is ion forward type?
  if(ion->type == B_ION ||
     ion->type == A_ION ||
     ion->type == C_ION)
    {
      return TRUE;
    }

  // reverse type ion
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

  // only add ions with no modifications
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
FLOAT_T get_ion_mass_z(
  ION_T* working_ion///< return the location of this ion -in 
  )
{
  if( working_ion == NULL ){
    return 0;
  }

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
  FLOAT_T mass_z ///< the m/z location -in
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

