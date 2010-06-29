/****************************************************************//**
 * \file alphabet.cpp
 * AUTHOR: Unknown
 * \brief Provide a function for converting amino acids to ints.
 ********************************************************************/
#include "alphabet.h"

static const int AMINO_ARRAY_CAPACITY = 100;
static int amino_array[AMINO_ARRAY_CAPACITY];

/**
 * Converts a character into an amino acid
 */
static void populate_amino_array(void){
  if (
    amino_array['T'-'A'] == 16 &&
    amino_array['N'-'A'] == 11
      ){
    return;
  }
  amino_array['A'-'A'] = 0;
  amino_array['C'-'A'] = 1;
  amino_array['D'-'A'] = 2;
  amino_array['E'-'A'] = 3;
  amino_array['F'-'A'] = 4;
  amino_array['G'-'A'] = 5;
  amino_array['H'-'A'] = 6;
  amino_array['I'-'A'] = 7;
  amino_array['K'-'A'] = 8;
  amino_array['L'-'A'] = 9;
  amino_array['M'-'A'] = 10;
  amino_array['N'-'A'] = 11;
  amino_array['P'-'A'] = 12;
  amino_array['Q'-'A'] = 13;
  amino_array['R'-'A'] = 14;
  amino_array['S'-'A'] = 15;
  amino_array['T'-'A'] = 16;
  amino_array['V'-'A'] = 17;
  amino_array['W'-'A'] = 18;
  amino_array['Y'-'A'] = 19;
  return;
}

/**
 * Converts a character into an int using the global amino_hash 
 */
int amino_to_int(char amino){
  populate_amino_array();
  int value = amino_array[amino - 'A'];
  return value;
}
