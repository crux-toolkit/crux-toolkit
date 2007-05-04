/**
 * \file mass.c 
 * $Revision: 1.8 $
 * \brief Provides constants and methods for calculating mass
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "objects.h"
#include "mass.h"
#include "protein.h"
#include "peptide.h"
#include "carp.h"
#include "ion.h"

/**
 * Array to store the amino acid masses
 */
float amino_masses[('Z' - 'A')*2 + 2];

/**
 * Have we initialized the amino acid masses?
 */
BOOLEAN_T initialized_amino_masses = FALSE;

//FIXME need to find the monoisotopic mass for some AA -chris
/**
 * initializes the mass array
 */
void initialize_amino_masses (void)
{
  //average mass
  amino_masses['A' - 'A'] = 71.0788;
  amino_masses['B' - 'A'] = 114.5962;
  amino_masses['C' - 'A'] = 103.1388;// + 57.000;
  amino_masses['D' - 'A'] = 115.0886;
  amino_masses['E' - 'A'] = 129.1155;
  amino_masses['F' - 'A'] = 147.1766;
  amino_masses['G' - 'A'] = 57.0519;
  amino_masses['H' - 'A'] = 137.1411;
  amino_masses['I' - 'A'] = 113.1594;
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
  amino_masses['X' - 'A'] = 113.1594;
  amino_masses['Y' - 'A'] = 163.1760;
  amino_masses['Z' - 'A'] = 128.6231;
  
  //monoisotopic mass
  amino_masses['A' - 'A' + 26] = 71.03711;
  amino_masses['B' - 'A' + 26] = 114.53494;
  amino_masses['C' - 'A' + 26] = 103.00919;// + 57.000;
  amino_masses['D' - 'A' + 26] = 115.02694;
  amino_masses['E' - 'A' + 26] = 129.04259;
  amino_masses['F' - 'A' + 26] = 147.06841;
  amino_masses['G' - 'A' + 26] = 57.02146;
  amino_masses['H' - 'A' + 26] = 137.05891;
  amino_masses['I' - 'A' + 26] = 113.08406;
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
  amino_masses['X' - 'A' + 26] = 113.08406;
  amino_masses['Y' - 'A' + 26] = 163.06333;
  amino_masses['Z' - 'A' + 26] = 128.55059;
  
}

/**
 * \returns The mass of the given amino acid.
 */
float get_mass_amino_acid(
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
    die("ERROR: mass type does not exist\n");
    //avoid compiler warning
    return -1;
  }
}

/**
 * \returns The average mass of the given amino acid.
 */
float get_mass_amino_acid_average(
  char amino_acid ///< the query amino acid -in
  )
{
  // has the amino_masses array been initialized?
  if(!initialized_amino_masses){
    initialize_amino_masses();
    initialized_amino_masses = TRUE;
  }
  return(amino_masses[(short int)amino_acid - 'A']);
}

/**
 * \returns The monoisotopic mass of the given amino acid.
 */
float get_mass_amino_acid_monoisotopic(
  char amino_acid ///< the query amino acid -in
  )
{
  // has the amino_masses array been initialized?
  if(!initialized_amino_masses){
    initialize_amino_masses();
    initialized_amino_masses = TRUE;
  }
  return(amino_masses[(short int)amino_acid - 'A' + 26 ]);
}

