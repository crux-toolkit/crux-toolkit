/**
 * \file mass.cpp 
 * \brief Provides constants and methods for calculating mass
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "objects.h"
#include "mass.h"
#include "Protein.h"
#include "Peptide.h"
#include "carp.h"
#include "Ion.h"
#include "modifications.h"

/* Private Variables */

/**
 * Array to store the amino acid masses
 */
FLOAT_T amino_masses[('Z' - 'A')*2 + 2];
enum {NUM_MOD_MASSES = 2048}; // 2 ^ MAX_AA_MODS = 2^11 = 2048
FLOAT_T aa_mod_masses[(int)NUM_MOD_MASSES];


/**
 * Have we initialized the amino acid masses?
 */
bool initialized_amino_masses = false;

// FIXME need to find the monoisotopic mass for some AA -chris
/**
 * initializes the mass array
 */
void initialize_amino_masses (void)
{
  // average mass
  amino_masses['A' - 'A'] = 71.0788;
  // set mass high to prevent peptides containing B
  amino_masses['B' - 'A'] = 7000.0;
  amino_masses['C' - 'A'] = 103.1388;
  amino_masses['D' - 'A'] = 115.0886;
  amino_masses['E' - 'A'] = 129.1155;
  amino_masses['F' - 'A'] = 147.1766;
  amino_masses['G' - 'A'] = 57.0519;
  amino_masses['H' - 'A'] = 137.1411;
  amino_masses['I' - 'A'] = 113.1594;
  amino_masses['J' - 'A'] = 113.1594;
  amino_masses['K' - 'A'] = 128.1741;
  amino_masses['L' - 'A'] = 113.1594;
  amino_masses['M' - 'A'] = 131.1926;
  amino_masses['N' - 'A'] = 114.1038;
  amino_masses['O' - 'A'] = 114.1472;
  amino_masses['P' - 'A'] = 97.1167;
  amino_masses['Q' - 'A'] = 128.1307;
  amino_masses['R' - 'A'] = 156.1875;
  amino_masses['S' - 'A'] = 87.0782;
  amino_masses['T' - 'A'] = 101.1051;
  amino_masses['U' - 'A'] = 150.0388;
  amino_masses['V' - 'A'] = 99.1326;
  amino_masses['W' - 'A'] = 186.2132;
  // set mass high to prevent peptides containing X
  amino_masses['X' - 'A'] = 7000.0;
  amino_masses['Y' - 'A'] = 163.1760;
  // set mass high to prevent peptides containing Z
  amino_masses['Z' - 'A'] = 7000.0;
  
  // monoisotopic mass
  amino_masses['A' - 'A' + 26] = 71.03711;
  // set mass high to prevent peptides containing B
  amino_masses['B' - 'A' + 26] = 7000.0;
  amino_masses['C' - 'A' + 26] = 103.00919;
  amino_masses['D' - 'A' + 26] = 115.02694;
  amino_masses['E' - 'A' + 26] = 129.04259;
  amino_masses['F' - 'A' + 26] = 147.06841;
  amino_masses['G' - 'A' + 26] = 57.02146;
  amino_masses['H' - 'A' + 26] = 137.05891;
  amino_masses['I' - 'A' + 26] = 113.08406;
  amino_masses['J' - 'A' + 26] = 113.08406;
  amino_masses['K' - 'A' + 26] = 128.09496;
  amino_masses['L' - 'A' + 26] = 113.08406;
  amino_masses['M' - 'A' + 26] = 131.04049;
  amino_masses['N' - 'A' + 26] = 114.04293;
  amino_masses['O' - 'A' + 26] = 114.07931;
  amino_masses['P' - 'A' + 26] = 97.05276;
  amino_masses['Q' - 'A' + 26] = 128.05858;
  amino_masses['R' - 'A' + 26] = 156.10111;
  amino_masses['S' - 'A' + 26] = 87.03203;
  amino_masses['T' - 'A' + 26] = 101.04768;
  amino_masses['U' - 'A' + 26] = 150.04344;
  amino_masses['V' - 'A' + 26] = 99.06841;
  amino_masses['W' - 'A' + 26] = 186.07931;
  // set mass high to prevent peptides containing X
  amino_masses['X' - 'A' + 26] = 7000.0;
  amino_masses['Y' - 'A' + 26] = 163.06333;
  // set mass high to prevent peptides containing Z
  amino_masses['Z' - 'A' + 26] = 7000.0;
  
  // modifications
  initialize_aa_mod_combinations_array();

  initialized_amino_masses = true;
}

/**
 * \brief Populates the array aa_mod_masses with the mass change of
 * all possible combinations of aa_mods.  Gets the list of aa_mods
 * from parameter.c.  For example, if mod1 and a mass change of 50 and
 * mod5 has a mass change of 10, then the entry at index 
 * (binary 00010001 =) 17 is (50 + 10=) 60.
 */
void initialize_aa_mod_combinations_array(){
  // set all to 0
  int i = 0;
  for(i = 0; i < NUM_MOD_MASSES; i++){
    aa_mod_masses[i] = 0;
  }

  // get aa mod list
  AA_MOD_T** amod_list = NULL;
  int num_mods = get_all_aa_mod_list(&amod_list);
  int mod_idx = 0;

  // for each aa mod
  for(mod_idx = 0; mod_idx < num_mods; mod_idx++){
    //printf("\nmod#: %i\n", mod_idx);
    int entry_idx = 0;

    while( entry_idx < NUM_MOD_MASSES ){
      //printf("skip from %i", entry_idx);
      // skip 2^i entries for mod_i (first mod is 1)
      entry_idx += (int)pow(2.0, mod_idx);
      //printf(" to %i\n", entry_idx);

      // write i entries for mod_i
      //int write_num = mod_idx+1;
      int write_num = (int)pow(2.0, mod_idx);
      while(write_num > 0 && entry_idx < NUM_MOD_MASSES){
        //printf("write at %i from %.3f", entry_idx, aa_mod_masses[entry_idx]);

        aa_mod_masses[entry_idx] += aa_mod_get_mass_change(amod_list[mod_idx]);
        //printf(" to  %.3f\n", aa_mod_masses[entry_idx]);
        entry_idx++;
        write_num--;
      }// next write
    }// next skip
  }// next mod

  //printf("Initialized mod combinations are:\n");
  //for(i=0; i<36; i++){
    //printf("index %i\tmass %.1f\n", i, aa_mod_masses[i]);
  //}
}

/**
 * \returns The mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid(
  char amino_acid, ///< the query amino acid -in
  MASS_TYPE_T mass_type ///< the isotopic mass type (AVERAGE, MONO) -in
  )
{
  if(mass_type == AVERAGE){
    return get_mass_amino_acid_average(amino_acid);
  }
  else if(mass_type == MONO){
    return get_mass_amino_acid_monoisotopic(amino_acid);
  }
  else{
    carp(CARP_FATAL, "ERROR: mass type does not exist\n");
    // avoid compiler warning
    return 1;
  }
}

/**
 * \returns The mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid(
  MODIFIED_AA_T amino_acid, ///< the query amino acid -in
  MASS_TYPE_T mass_type ///< the isotopic mass type (AVERAGE, MONO) -in
  ){

  switch( mass_type ){
  case AVERAGE:
    return get_mass_mod_amino_acid_average(amino_acid);
  case MONO:
    return get_mass_mod_amino_acid_monoisotopic(amino_acid);
  case NUMBER_MASS_TYPES:
    carp(CARP_FATAL,"Bad mass type.");
  }
  return 0;
}

/**
 * \returns The average mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid_average(
  char amino_acid ///< the query amino acid -in
  )
{
  // has the amino_masses array been initialized?
  if(!initialized_amino_masses){
    initialize_amino_masses();
    initialized_amino_masses = true;
  }
  return(amino_masses[(short int)amino_acid - 'A']);
}

/**
 * \returns The average mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid_average(
  MODIFIED_AA_T amino_acid ///< the query amino acid -in
  ){
  if(!initialized_amino_masses){
    initialize_amino_masses();
  }
  //printf("aa is %hu, char %c\n", (unsigned short)amino_acid, modified_aa_to_char(amino_acid));

  short int aa = amino_acid & GET_AA_MASK;
  unsigned short int mod = amino_acid & GET_MOD_MASK;
  mod = mod >> 5;
  //printf("aa idx %hi mass %.3f + mod idx %hu mass %.3f\n", aa,  amino_masses[aa], mod, aa_mod_masses[mod]);
  return amino_masses[aa] + aa_mod_masses[mod];

}

/**
 * Finds the modification identifier associated with the given mass
 * shift.  Can be the identifier from a single modfification or from
 * multiple modficiations to the same residue.  The returned
 * identifier can be used to modify a MODIFIED_AA_T so that it has the
 * given mass shift. 
 */
MODIFIED_AA_T get_mod_identifier(FLOAT_T mass_shift){

  if(!initialized_amino_masses){
    initialize_amino_masses();
  }

  int precision = get_int_parameter("mod-precision");
  for(int mod_idx = 0; mod_idx < (int)NUM_MOD_MASSES; mod_idx++){
    if( is_equal(mass_shift, aa_mod_masses[mod_idx], precision) ){
      MODIFIED_AA_T identifier = mod_idx;
      identifier = identifier << 5;
      return identifier;
    }
  }

  // if we got to here, no mod found
  carp_once(CARP_WARNING, "No modification identifier found for mass shift %.*f."
                          "\nWarning Suppressed, others may exist",
       precision, mass_shift);

  return 0;
}

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
FLOAT_T get_mass_amino_acid_monoisotopic(
  char amino_acid ///< the query amino acid -in
  )
{
  // has the amino_masses array been initialized?
  if(!initialized_amino_masses){
    initialize_amino_masses();
    initialized_amino_masses = true;
  }
  return(amino_masses[(short int)amino_acid - 'A' + 26 ]);
}

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
FLOAT_T get_mass_mod_amino_acid_monoisotopic(
  MODIFIED_AA_T amino_acid ///< the query amino acid -in
  ){
  if(!initialized_amino_masses){
    initialize_amino_masses();
  }

  short int aa = amino_acid & GET_AA_MASK;
  unsigned short int mod = amino_acid & GET_MOD_MASK;
  mod = mod >> 5;
  //printf("mod aa idx %hi mass %.3f + mod idx %hu mass %.3f\n", aa,  amino_masses[aa+26], mod, aa_mod_masses[mod]);
  return amino_masses[aa + 26] + aa_mod_masses[mod];

}


/**
 * increase the amino acid mass for both mono and average
 */
void increase_amino_acid_mass(
  char amino_acid, ///< the query amino acid -in
  FLOAT_T update_mass ///< the mass amount to update for the amino acid -in
  )
{
  // has the amino_masses array been initialized?
  if(!initialized_amino_masses){
    initialize_amino_masses();
    initialized_amino_masses = true;
  }

  // check if amino acid
  if((short int)amino_acid < 'A' || (short int)amino_acid > 'Z'){
    carp(CARP_ERROR, "Cannot update mass, char: %c not an amino acid");
  }
  
  amino_masses[(short int)amino_acid - 'A' + 26 ] += update_mass;
  amino_masses[(short int)amino_acid - 'A'] += update_mass;
}
