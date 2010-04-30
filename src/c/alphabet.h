/**
 * \file alphabet.h
 * AUTHOR: Unknown
 * \brief Provide a function for converting amino acids to ints.
 ********************************************************************/
#ifndef ALPHABET_H
#define ALPHABET_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Converts a character into an int using the global amino_hash 
 */
int amino_to_int
  (char amino);

#ifdef __cplusplus
}
#endif

#endif
