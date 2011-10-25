/**
 * \file alphabet.h
 * AUTHOR: Unknown
 * \brief Provide a function for converting amino acids to ints.
 ********************************************************************/
#ifndef ALPHABET_H
#define ALPHABET_H


class Alphabet {

 public:
  /**
   * Converts a character into an int using the global amino_hash 
   */
  static int aminoToInt(char amino);

  /**
   * Constants for iterating through all amino acids.
   */
  static const char* aminoAcids[];
  static const int numAminoAcids;

 protected:
  static const int AMINO_ARRAY_CAPACITY = 100;
  static int amino_array_[AMINO_ARRAY_CAPACITY];
  static bool amino_array_populated_;
  static void populateAminoArray();

};
#endif
