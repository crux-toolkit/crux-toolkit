/****************************************************************//**
 * \file alphabet.cpp
 * AUTHOR: Unknown
 * \brief Provide a function for converting amino acids to ints.
 ********************************************************************/
#include "Alphabet.h"

int Alphabet::amino_array_[AMINO_ARRAY_CAPACITY] = {0};
bool Alphabet::amino_array_populated_ = false;

/**
 * Converts a character into an amino acid
 */
void Alphabet::populateAminoArray(void){
  if ( amino_array_populated_ ){
    return;
  }
  amino_array_['A'-'A'] = 0;
  amino_array_['C'-'A'] = 1;
  amino_array_['D'-'A'] = 2;
  amino_array_['E'-'A'] = 3;
  amino_array_['F'-'A'] = 4;
  amino_array_['G'-'A'] = 5;
  amino_array_['H'-'A'] = 6;
  amino_array_['I'-'A'] = 7;
  amino_array_['K'-'A'] = 8;
  amino_array_['L'-'A'] = 9;
  amino_array_['M'-'A'] = 10;
  amino_array_['N'-'A'] = 11;
  amino_array_['P'-'A'] = 12;
  amino_array_['Q'-'A'] = 13;
  amino_array_['R'-'A'] = 14;
  amino_array_['S'-'A'] = 15;
  amino_array_['T'-'A'] = 16;
  amino_array_['V'-'A'] = 17;
  amino_array_['W'-'A'] = 18;
  amino_array_['Y'-'A'] = 19;

  amino_array_populated_ = true;
  return;
}

/**
 * Converts a character into an int using the global amino_hash 
 */
int Alphabet::aminoToInt(char amino){
  populateAminoArray();
  int value = amino_array_[amino - 'A'];
  return value;
}


const char* Alphabet::aminoAcids[]  = { "A", "C", "D", "E", "F", "G", "H", 
                                        "I", "K", "L", "M", "N", "P", "Q", 
                                        "R", "S", "T", "V", "W", "Y", };
const int Alphabet::numAminoAcids = sizeof(aminoAcids) / sizeof(char*);
