/*************************************************************************//**
 * \file ion.cpp
 * \brief Object for representing a single ion.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "objects.h"
#include "Ion.h"
#include "Alphabet.h"
#include "Peptide.h"
#include "Peak.h"
#include "mass.h"
#include "utils.h"
#include <sys/types.h>
#ifndef _MSC_VER
#include <netinet/in.h>
#include <inttypes.h>
#else
#include <Winsock2.h>
#include "WinCrux.h"
#endif

using namespace Crux;

// At one point I need to reverse the endianness for pfile_create to work
  // Apparently that is no longer true. Hence 0 below.
  static const bool REVERSE_ENDIAN = 0;

  /**
   * Array to store the modification masses
   */
  static FLOAT_T modification_masses[MAX_MODIFICATIONS];

  /**
   * Have we initialized the modification_masses?
   */
  static bool initialized_modification_masses = false;

  static const int DETECTABLE_MZ_MIN = 200;
  static const int DETECTABLE_MZ_MAX = 2400;


/**
 * initializes the mass array
 */
void Ion::initializeModificationMasses(
   MASS_TYPE_T mass_type ///< mass type (average, mono) -in
)
{
  // set modification mass
  if(mass_type == MONO){
    modification_masses[NH3] = -MASS_NH3_MONO;
    modification_masses[H2O] = -MASS_H2O_MONO ;
    modification_masses[ISOTOPE] = 1; // FIXME check this!!!
    modification_masses[FLANK] = 1; // FIXME check this!!!
  }
  else if(mass_type == AVERAGE){
    modification_masses[NH3] = -MASS_NH3_AVERAGE;
    modification_masses[H2O] = -MASS_H2O_AVERAGE;
    modification_masses[ISOTOPE] = 1; // FIXME check this!!!
    modification_masses[FLANK] = 1; // FIXME check this!!!
  }

  initialized_modification_masses = true;
}

/**
 * initializes an Ion object.
 */
void Ion::init() {
  type_ = (ION_TYPE_T)0;
  cleavage_idx_ = 0;
  charge_ = 0;
  peptide_sequence_ = NULL;
  peptide_mass_ = 0;
  memset(modification_counts_,0,sizeof(int)*MAX_MODIFICATIONS);
  ion_mass_z_ = 0;
  peak_ = NULL;
  pointer_count_ = 0;
}

/**
 * \returns An (empty) ion object.
 */
Ion::Ion() {
  init();
}

/**
 * helper function
 * only copies the pointer to the peptide sequence
 * creates an heap allocated ion, ion mass is not calculated.
 * \returns an ION_T object
 */
void Ion::initBasicIon(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide ///< location for the new ion -in
  ) 
{
  // init ion
  init();
  type_ = type;
  cleavage_idx_ = cleavage_idx;
  charge_ = charge;
  peptide_sequence_ = peptide;
  // TODO get mass type from param file
  peptide_mass_ = Peptide::calcSequenceMass(peptide, MONO); 
  peak_ = NULL;
}

/**
 * only copies the pointer to the peptide sequence
 * ion mass/z with out any modification
 * \returns an ION_T object
 */
Ion::Ion(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // get new basic ion
  initBasicIon(type, cleavage_idx, charge, peptide);
  
  // calculate and set ion mass/z
  if(!calcMassZ(mass_type, false)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
  }
}

/**
 * only copies the pointer to the peptide sequence
 * inputs a array of all the modification counts
 * \returns an ION_T object
 */
Ion::Ion(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  int* modification_counts ///< an array of modification counts for each modification -in
  )
{
  // get new basic ion
  initBasicIon(type, cleavage_idx, charge, peptide);
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    modification_counts_[modification_idx] = modification_counts[modification_idx];
  }
  
  // calculate and set ion mass/z
  if(!calcMassZ(mass_type, true)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
  }
}

/**
 * only copies the pointer to the peptide sequence
 * inputs a array of all the modification counts
 * inputs the pre modified mass, of just all AA mass summed up.
 * \returns an ION_T object
 */
Ion::Ion(
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
  initBasicIon(type, cleavage_idx, charge, peptide);
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    modification_counts_[modification_idx] = modification_counts[modification_idx];
  }
  
  // calculate and set ion mass/z
  if(!calcMassZWithMass(mass_type, base_mass, true)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
  } 
}

/**
 * only copies the pointer to the peptide sequence
 * inputs the pre modified mass, of just all AA mass summed up.
 * \returns an ION_T object
 */
Ion::Ion(
  ION_TYPE_T type,   ///< intensity for the new ion -in 
  int cleavage_idx, ///< index into the peptide amide bonds of this ion
  int charge, ///< charge of the ion
  char* peptide, ///< location for the new ion -in
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  FLOAT_T base_mass ///< the base mass of the ion -in
  )
{
  // get new basic ion
  initBasicIon(type, cleavage_idx, charge, peptide);

  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    modification_counts_[modification_idx] = 0;
  }
  
  // calculate and set ion mass/z
  if(!calcMassZWithMass(mass_type, base_mass, true)){
    carp(CARP_ERROR, "failed to calculate ion mass/z");
  } 
}


/**
 * frees A ION_T object
 */
Ion::~Ion() {
}

/**
 * decrements the pointer and frees ion 
 */
void Ion::freeIon(
  Ion* ion ///< the ion to free
  ) {

  ion->pointer_count_--;

  if (ion->pointer_count_ <= 0) {
    delete ion;
  }
}

/**
 * increments the pointer count for the ion
 */
void Ion::incrementPointerCount() {
  pointer_count_++;
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
void Ion::print(
  FILE* file ///< to this file -in
  )
{
  // print all fields of ion


  char* type_str = ion_type_to_string(type_);

  fprintf(file, "%.2f\t%.2f\t%d\t%s\t%d", 
      ion_mass_z_, (ion_mass_z_)*charge_, charge_, 
      type_str, cleavage_idx_);

  free(type_str);

  // iterate over all modification counts
  int mod_idx;
  for(mod_idx=0; mod_idx < MAX_MODIFICATIONS; ++mod_idx){
    fprintf(file,"\t%d", modification_counts_[mod_idx]);
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
void Ion::printGmtkSingle(
  FILE* file  ///< to this file -in
  ){

  int has_mobile_proton = 1;
  int is_detectable = 0;
  int is_detected = 0;
  int is_possible = 0;

  FLOAT_T intensity = 0.0;
  FLOAT_T intensity_rank = 0.0;
  if (peak_ != NULL){
    intensity = peak_->getIntensity();
    intensity_rank = peak_->getIntensityRank();
    is_detected = 1;
  }

  if ((ion_mass_z_ >= DETECTABLE_MZ_MIN) &&
      (ion_mass_z_ <= DETECTABLE_MZ_MAX)){
    is_detectable = 1;
  }

  FLOAT_T mz_ratio = (ion_mass_z_)/(peptide_mass_);
  int mz_int = (int)(mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);

  const char* format = "%.6f\t%.6f\t%.6f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n";
  fprintf(file, format,
      mz_ratio,                                                 // 1 
      intensity,                                                // 2 
      intensity_rank,                                           // 3 
      has_mobile_proton,                                        // 4 
      mz_int,                                                   // 5 
      cleavage_idx_,                                            // 6
      strlen(peptide_sequence_) - cleavage_idx_ + 1,            // 7
      Alphabet::aminoToInt(peptide_sequence_[cleavage_idx_-1]),// 8 
      Alphabet::aminoToInt(peptide_sequence_[cleavage_idx_]), // 9 
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
void Ion::printGmtkSingleBinary(
  FILE* file,  ///< to this file -in
  int sentence_idx,
  int frame_idx
  ){

  FLOAT_T* float_array = (FLOAT_T*)mycalloc(sizeof(FLOAT_T), SINGLE_ION_FLOATS);
  int* int_array = (int*)mycalloc(sizeof(int), SINGLE_ION_INTS);

  FLOAT_T mz_ratio = (ion_mass_z_)/(peptide_mass_);
  float_array[0] = mz_ratio;                              // 0
  float_array[1] = 0.0;                                   // 1
  float_array[2] = 0.0;                                   // 2
  int is_detected = 0;
  if (peak_ != NULL){
    float_array[1] = peak_->getIntensity();      //1
    float_array[2] = peak_->getIntensityRank();  //2
    is_detected = 1; 

#ifdef LOG_INTENSITY
    // START
    // replace the rank intensity with the negative log of the 
    // TIC normalized intensity. A hack, I know
    float_array[2] = 10.0 + log(float_array[1]);
#endif

  }

  int mz_int = (int)(mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int cterm_idx = strlen(peptide_sequence_) - cleavage_idx_; 
  int left_amino = Alphabet::aminoToInt(peptide_sequence_[cleavage_idx_-1]);
  int right_amino = Alphabet::aminoToInt(peptide_sequence_[cleavage_idx_]);
  int is_detectable = 0;
  if ( 
       ((ion_mass_z_ >= DETECTABLE_MZ_MIN) 
              &&
        (ion_mass_z_ <= DETECTABLE_MZ_MAX)) 

         || 
         
       is_detected

     ){
    is_detectable = 1;                                    
  }


  int_array[0] = 1;                 // 0
  int_array[1] = mz_int;            // 1
  int_array[2] = cleavage_idx_;     // 2
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
void Ion::printNullGmtkSingleBinary(
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
void Ion::printNullGmtkPairedBinary(
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
void Ion::printNullGmtkSingle(
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
void Ion::printGmtkPairedBinary(
  Ion* first_ion,   ///< print this ion -in
  Ion* second_ion,  ///< print this ion -in
  FILE* file,         ///< to this file -in
  int sentence_idx, 
  int frame_idx
  ){
  
  FLOAT_T* float_array = (FLOAT_T*)mycalloc(sizeof(FLOAT_T), PAIRED_ION_FLOATS);
  int* int_array = (int*)mycalloc(sizeof(int), PAIRED_ION_INTS);

  // start with the floats
  FLOAT_T n_mz_ratio = (first_ion->ion_mass_z_)/(first_ion->peptide_mass_);
  float_array[0] = n_mz_ratio;                                    // 0
  // TODO 
  // subtract from 1.0?
  FLOAT_T c_mz_ratio =  1.0 - n_mz_ratio;                           // 1
  float_array[1] = c_mz_ratio;

  int first_is_detected = 0;
  if (first_ion->peak_ != NULL){
    // put in LOG_INTENSITY ?
    float_array[2] = first_ion->peak_->getIntensity();    //2
    float_array[4] = first_ion->peak_->getIntensityRank();//4
    first_is_detected = 1; 
  }

  int second_is_detected = 0;
  if (second_ion->peak_ != NULL){
    // put in LOG_INTENSITY ?
    float_array[3] = second_ion->peak_->getIntensity();    //3
    float_array[5] = second_ion->peak_->getIntensityRank();//5
    second_is_detected = 1; 
  }

  // next do the ints
  int n_mz_int = (int)(n_mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int c_mz_int = (int)(c_mz_ratio * (MZ_INT_MAX - MZ_INT_MIN) + MZ_INT_MIN);
  int cterm_idx = strlen(first_ion->peptide_sequence_) 
    - first_ion->cleavage_idx_; 
  int left_amino = Alphabet::aminoToInt(
      first_ion->peptide_sequence_[first_ion->cleavage_idx_-1]);
  int right_amino = Alphabet::aminoToInt(
      first_ion->peptide_sequence_[first_ion->cleavage_idx_]);
  int first_is_detectable = 0;
  if ( 
       ((first_ion->ion_mass_z_ >= DETECTABLE_MZ_MIN) 
              &&
        (first_ion->ion_mass_z_ <= DETECTABLE_MZ_MAX)) 

         || 
         
       first_is_detected

     ){
    first_is_detectable = 1;                                    
  }

  int second_is_detectable = 0;
  if ( 
       ((second_ion->ion_mass_z_ >= DETECTABLE_MZ_MIN) 
              &&
        (second_ion->ion_mass_z_ <= DETECTABLE_MZ_MAX)) 

         || 
         
       second_is_detected

     ){
    second_is_detectable = 1;                                    
  }

  int_array[0] = n_mz_int; 
  int_array[1] = c_mz_int; 
  int_array[2] = first_ion->cleavage_idx_;
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
FLOAT_T Ion::modifyMassZ(
  FLOAT_T mass_z, ///< the pre-modified ion mass_z -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  int charge, ///< charge of the ion
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // initalize the modified masses(average|mono);
  if(!initialized_modification_masses){
    initializeModificationMasses(mass_type);
  }
  
  return  mass_z + (modification_masses[(int)ion_modification]/(float)charge) * modification_count;  
}


/**
 * Adds the given ION_MODIFICATION to this ion
 */
void Ion::addModification(
  ION_MODIFICATION_T modification, ///< add this modification to the ion -in
  int modification_count, ///< the number of modifications -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // update modification count
  modification_counts_[(int)modification] += modification_count;  
  // reset ion mass_z
  ion_mass_z_ =  modifyMassZ(ion_mass_z_, 
                                   modification_count, 
                                   modification, 
                                   charge_, 
                                   mass_type);
}

/**
 * Adds the given ION_MODIFICATION to this ion
 */
void Ion::setPeak(
  Peak * peak ///< peak to add to this ion -in
  ){
  peak_ = peak;
};

/**
 *\returns the ion's AA mass added all up
 */
FLOAT_T Ion::getMass(
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )                   
{
  char* ion_sequence = NULL;
  ION_TYPE_T ion_type = type_;
  int ion_length = 0;
  bool memory_used = false;
  FLOAT_T mass = 0;
  bool reverse = false;

  // get sequence for x,y,z ion
  if(ion_type == X_ION ||ion_type == Y_ION || ion_type == Z_ION){
    // convert the cleavage index into the actually index that start from the left.
    int real_cleavage_idx = strlen(peptide_sequence_) - cleavage_idx_;
    ion_sequence = &(peptide_sequence_[real_cleavage_idx]);
    ion_length = strlen(ion_sequence);
    reverse = true;
  }
  // get sequence for a,b,c ion
  else{
    ion_length = cleavage_idx_;
    ion_sequence = (char*)mycalloc(ion_length+1, sizeof(char));
    strncpy(ion_sequence, peptide_sequence_, ion_length);
    memory_used = true;
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
    }else{
      mass += MASS_H2O_MONO;
    }
  }
  
  return mass;
}

/**
 *\return the modified mass accodring to the modification type
 */
FLOAT_T Ion::modifyMass(
  FLOAT_T mass, ///< the pre-modified ion mass -in
  int modification_count, ///< number of times more the modification occurs -in
  ION_MODIFICATION_T ion_modification, ///< the type of modification -in
  MASS_TYPE_T mass_type ///< mass type (average, mono) -in
  )
{
  // initalize the modified masses(average|mono);
  if(!initialized_modification_masses){
    initializeModificationMasses(mass_type);
  }
  
  return  mass + modification_masses[(int)ion_modification] * modification_count;
}

/**
 * is_modified, indiciates if there are any modification to the ion
 * speeds up the proccess if FALSE.
 *\returns true if successfully computes the mass/z of the ion, else FALSE
 */
bool Ion::calcMassZWithMass(
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  FLOAT_T mass, ///< the basic mass of the ion -in
  bool is_modified ///< are there any modifications for this ion? -in
  )
{
  FLOAT_T h_mass = MASS_H_MONO;

  // alter mass according to the modification
  if(is_modified){
    // iterate over all type of modifications for ion
    int modification_idx = 0;
    for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
      // update ion mass if modification is needed
      if(modification_counts_[modification_idx] != 0){
        mass = modifyMass(mass, modification_counts_[modification_idx], 
                               (ION_MODIFICATION_T)modification_idx, mass_type);
      }
    }
  }
  
  // reset h mass if needed
  if(mass_type == AVERAGE){
    h_mass = MASS_H_AVERAGE;
  }

  // convert mass to m/z and assigned to ion
  ion_mass_z_ = (mass + (h_mass*(FLOAT_T)charge_))/(FLOAT_T)charge_;

  return true;
}

/**
 * is_modified, indiciates if there are any modification to the ion
 * speeds up the proccess if FLASE.
 *\returns true if successfully computes the mass/z of the ion, else FALSE
 */
bool Ion::calcMassZ(
  MASS_TYPE_T mass_type, ///< mass type (average, mono) -in
  bool is_modified ///< are there any modifications for this ion? -in
  )
{
  FLOAT_T mass = this->getMass(mass_type);
  return this->calcMassZWithMass(mass_type, mass, is_modified);
}

/**
 * Copies ion object from src to dest.
 * must pass in a seperate pointer peptide sequence from its own ion_series object
 * must pass in a memory allocated ION_T* dest
 */
void Ion::copy(
  Ion* src,///< ion to copy from -in
  Ion* dest,///< ion to copy to -out
  char* peptide_sequence ///< the peptide sequence that the dest should refer to -in
  )
{
  dest->type_ = src->type_;
  dest->cleavage_idx_ = src->cleavage_idx_;
  dest->charge_ = src->charge_;
  
  // set all modification counts in the ion
  int modification_idx = 0;
  for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
    dest->modification_counts_[modification_idx] = src->modification_counts_[modification_idx];
  }
  
  dest->ion_mass_z_ = src->ion_mass_z_;
  dest->peptide_sequence_ = peptide_sequence;
  dest->peptide_mass_ = src->peptide_mass_;
}


/**
 * \returns true if forward ion_type(A,B,C), 
 * else reverse ion_type(X,Y,Z) FALSE
 */
bool Ion::isForwardType()
{
  // is ion forward type?
  if(type_ == B_ION ||
     type_ == A_ION ||
     type_ == C_ION)
    {
      return true;
    }

  // reverse type ion
  return false;
}

/**
 *\returns true if the ion has modifications, else FALSE
 */
bool Ion::isModified()
{
  int by_modification = 0;

  // only add ions with no modifications
  for(; by_modification < MAX_MODIFICATIONS; ++by_modification){
    if(modification_counts_[by_modification] != 0){
      return true;
    }
  }
  
  return false;
}

/*********************************
 * get, set methods for ion fields
 *********************************/

/**
 * \returns the location of ION_T object
 */
FLOAT_T Ion::getMassZ()
{
  return ion_mass_z_;
}

/**
 * sets the mass/z of the ION_T object
 * while this can be calculated from the char*, cleavage_idx and
 * modifications, it allows some optimizations if we allow it to be set
 * instead
 */
void Ion::setMassZ(
  FLOAT_T mass_z ///< the m/z location -in
  )
{
  ion_mass_z_ = mass_z;
}

FLOAT_T Ion::getMassFromMassZ() {
  FLOAT_T ans = (ion_mass_z_ - MASS_PROTON) * (FLOAT_T)charge_;
  return ans;
  
}

void Ion::setMassZFromMass(
  FLOAT_T mass
  ) {

  FLOAT_T charge = charge_;
  ion_mass_z_ = (mass + MASS_PROTON * charge) / charge;

}


/**
 * return the cleavage_idx of the ion object
 */
int Ion::getCleavageIdx()
{
  return cleavage_idx_;
}

/**
 * set the cleavage_idx of the ion object
 */
void Ion::setCleavageIdx(
  int cleavage_idx ///< the cleavage index in the peptide -in
  )
{
  cleavage_idx_ = cleavage_idx;
}

/**
 * return the charge of the ion object
 */
int Ion::getCharge()
{
  return charge_;
}

/**
 * set the charge of the ion object
 */
void Ion::setCharge(
  int charge ///< the charge of this ion -in
  )
{
  charge_ = charge;
}

/**
 * return the ION_TYPE_T of the ion object
 */
ION_TYPE_T Ion::getType()
{
  return type_;
}

/**
 * set the ION_TYPE_T of the ion object
 */
void Ion::setType(
  ION_TYPE_T ion_type ///< the ion type of this ion -in 
  )
{
  type_ = ion_type;
}

/**
 * return the parent peptide sequence of the ion object
 * returns a pointer to the sequence, should not free
 */
char* Ion::getPeptideSequence()
{
  return peptide_sequence_;
}

/**
 * return a pointer to the modification_count array of the ion object
 */
int* Ion::getModificationCounts()
{
  return modification_counts_;
}

/**
 * return the count of in the modification_count array of the ion object
 */
int Ion::getSingleModificationCount(
  ION_MODIFICATION_T mod_type ///< the modification count wanted -in
  )
{
  return modification_counts_[mod_type];
}

int Ion::getTotalModificationCount() {
  int total = 0;

  for (int idx=NH3;idx < ALL_MODIFICATION;idx++) {
    total += modification_counts_[idx];
  }

  return total;

}


/**
 * set the parent peptide_sequence of the ion object
 */
void Ion::setPeptideSequence(
  char* peptide_sequence ///< the parent peptide's sequence of this ion -in 
  )
{
  peptide_sequence_ = peptide_sequence;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

