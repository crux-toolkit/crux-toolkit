/****************************************************************//**
 * \file alphabet.cpp
 * FILE: alphabet.cpp
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4-17-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2001, Columbia University
 * DESCRIPTION: Define the amino acid and nucleotide alphabets.
 ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "alphabet.h"
#include "utils.h"
#include "carp.h"
#include "hash.h"

#ifndef AMBIG_DEBUG
#define AMBIG_DEBUG 0
#endif

#define AMINO_ARRAY_CAPACITY 100

/* Define global variables to hold alphabet info. */
ALPH_T alph = INVALID_ALPH;
int alph_size = 0;
int ambigs = 0;
char alphabet[NUM_AMINOS + NUM_AMINO_AMBIGS + 1];
char any_char = '\0';

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void set_alphabet
  (VERBOSE_T verbose,
   const char*     given_alphabet)
{
  /* Is the alphabet already initialized? */
  if (which_alphabet() != INVALID_ALPH) {
    assert((unsigned)alph_size == strlen(given_alphabet));
  } else if (strlen(given_alphabet) == NUM_AMINOS ||
             strlen(given_alphabet) == NUM_AMINOS+1) {
    alph = PROTEIN_ALPH;
    alph_size = NUM_AMINOS;
    ambigs = NUM_AMINO_AMBIGS;
    any_char = ANY_AMINO;
    strcpy(alphabet, given_alphabet);
    strcat(alphabet, AMINO_AMBIGS);
    /* chris edited
    if (verbose >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using amino acid alphabet (%s).\n", alphabet);
    }
    */
    carp(CARP_INFO,"Using amino acid alphabet (%s).", alphabet);
  }

  else if (strlen(given_alphabet) == NUM_BASES) {
    alph = DNA_ALPH;
    alph_size = NUM_BASES;
    ambigs = NUM_BASE_AMBIGS;
    any_char = ANY_BASE;
    strcpy(alphabet, given_alphabet);
    strcat(alphabet, BASE_AMBIGS);
    if (verbose >= NORMAL_VERBOSE) {
      carp(CARP_INFO, "Using nucleotide alphabet (%s).\n", alphabet);
    }
  } else {
    carp(CARP_FATAL, "Unrecognized alphabet (%s).\n", given_alphabet);
  }
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
ALPH_T which_alphabet
  (void)
{
  return(alph);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
int get_alph_size
  (const ALPH_SIZE_T which_size)
{
  if (alph == INVALID_ALPH) {
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  }

  switch (which_size) {
  case ALPH_SIZE :
    return(alph_size);
    break;
  case ALL_SIZE :
    return(alph_size + ambigs);
    break;
  case AMBIG_SIZE :
    return(ambigs);
  default :
    carp(CARP_FATAL, "Illegal alphabet size request.\n");
  }
  /* Unreachable. */
  return(alph_size);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
char * get_alphabet
  (BOOLEAN_T include_ambigs)
{
  static char short_alphabet[NUM_AMINOS + 1];
  static BOOLEAN_T first_time = TRUE;

  if (alph == INVALID_ALPH) {
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  }

  if (!include_ambigs) {
    if (first_time) {
      strncpy(short_alphabet, alphabet, alph_size);
      short_alphabet[alph_size] = '\0';
      first_time = FALSE;
    }
    return(short_alphabet);
  }
  return(alphabet);
}

/********************************************************************
 * See .h file for description
 ********************************************************************/
char get_alph_char
  (const int char_index)
{
  if (alph == INVALID_ALPH) {
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  } 
  if ((char_index < 0) || (char_index > get_alph_size(ALPH_SIZE))) {
    carp(CARP_FATAL, "Requested character outside of alphabet (%d).\n", char_index);
  }
  return(alphabet[char_index]);
}

char get_any_char
  (void)
{
  switch(alph) {
  case INVALID_ALPH :
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
    break;
  case PROTEIN_ALPH :
    return(ANY_AMINO);
    break;
  case DNA_ALPH :
    return(ANY_BASE);
    break;
  }
  /* Unreachable. */
  return(ANY_BASE);
}


/********************************************************************
 * See .h file for description
 ********************************************************************/
#ifndef ALPH_DEBUG
#define ALPH_DEBUG 0
#endif
int alphabet_index
  (const char   letter,    /* The letter we want to find an index for. */
   const char * alphabet)  /* Alphabet in which we're looking for an index. */ 
{
  static BOOLEAN_T first_time = TRUE; /* Is this the first call? */
  static int alph_indices[27];        /* Inverse alphabet array. */
  int        alph_size;               /* The length of the given alphabet. */
  int        i_alph;                  /* Alphabet index. */
  int        return_value;

  /* The first time the function is called, compute an inverse array. */
  if (first_time) {
    first_time = FALSE;

    /* First fill the array with -1s. */
    for (i_alph = 0; i_alph < 26; i_alph++)
      alph_indices[i_alph] = -1;

    /* Compute the position of each letter in the alphabet. */
    alph_size = strlen(alphabet);
    for (i_alph = 0; i_alph < alph_size; i_alph++)
      alph_indices[(int)(alphabet[i_alph] - 'A')] = i_alph;

    DEBUG_CODE
      (ALPH_DEBUG,
       fprintf(stderr, "alphabet=%s\n", alphabet);
       for (i_alph = 0; i_alph < 26; i_alph++)
       fprintf(stderr, "%d ", alph_indices[i_alph]);
       fprintf(stderr, "\n");
       );
  }


  /* Retrieve the index from the inverse array. */
  return_value = alph_indices[(int)(letter - 'A')];

  //  Die if the character is not a letter.
  if (return_value == -1) {
    carp(CARP_FATAL, "Non-alphabetic character (%c).\n", letter);
  }

  return(return_value);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void get_nrdb_frequencies
  (ARRAY_T* freqs) 
{
  if (which_alphabet() == PROTEIN_ALPH) {
    set_array_item( 0, 0.073164, freqs); /* A */
    set_array_item( 1, 0.018163, freqs); /* C */
    set_array_item( 2, 0.051739, freqs); /* D */
    set_array_item( 3, 0.062340, freqs); /* E */
    set_array_item( 4, 0.040283, freqs); /* F */
    set_array_item( 5, 0.069328, freqs); /* G */
    set_array_item( 6, 0.022428, freqs); /* H */
    set_array_item( 7, 0.056282, freqs); /* I */
    set_array_item( 8, 0.058493, freqs); /* K */
    set_array_item( 9, 0.091712, freqs); /* L */
    set_array_item(10, 0.023067, freqs); /* M */
    set_array_item(11, 0.046077, freqs); /* N */
    set_array_item(12, 0.050674, freqs); /* P */
    set_array_item(13, 0.040755, freqs); /* Q */
    set_array_item(14, 0.051897, freqs); /* R */
    set_array_item(15, 0.073802, freqs); /* S */
    set_array_item(16, 0.059411, freqs); /* T */
    set_array_item(17, 0.064362, freqs); /* V */
    set_array_item(18, 0.013341, freqs); /* W */
    set_array_item(19, 0.032682, freqs); /* Y */
  } else if (which_alphabet() == DNA_ALPH) {
    set_array_item( 0, 0.281774, freqs); /* A */
    set_array_item( 1, 0.222020, freqs); /* C */
    set_array_item( 2, 0.228876, freqs); /* G */
    set_array_item( 3, 0.267330, freqs); /* T */
  } else {
    carp(CARP_FATAL, "Illegal alphabet size (%d).\n", alph_size);
  }
  fill_in_ambiguous_chars(FALSE, freqs);
}


/********************************************************************
 * Make one position in an alphabetic frequency array the average of
 * a set of other positions.
 ********************************************************************/
static void set_ambiguity
  (BOOLEAN_T log_space,
   const char      target,
   const char*     sources,
   ARRAY_T*  freqs)
{
  int    i_source;
  int    num_sources;
  char   source;
  PROB_T source_count;
  PROB_T this_count;

  /* Add up the frequencies of the source characters. */
  source_count = 0.0;
  num_sources = strlen(sources);
  for (i_source = 0; i_source < num_sources; i_source++) {
    source = sources[i_source];
    this_count = get_array_item(alphabet_index(source, alphabet), freqs);

    if (log_space) {
      source_count = LOG_SUM(source_count, this_count);
    } else {
      source_count += this_count;
    }
  }

  /* Divide by the number of source characters. */
  if (log_space) {
    source_count -= my_log2(num_sources);
  } else {
    source_count /= num_sources;
  }

  /* Store the result. */
  set_array_item(alphabet_index(target, alphabet), source_count, freqs);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void fill_in_ambiguous_chars
  (BOOLEAN_T log_space,
   ARRAY_T* freqs)  /* The emission distribution to be extended. 
                       (Must be pre-mallocked large enough to accept
                       ambiguous characters). */
{
  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    set_ambiguity(log_space, 'B', "DN", freqs);
    set_ambiguity(log_space, 'Z', "EQ", freqs);
    set_ambiguity(log_space, 'U', "ACDEFGHIKLMNPQRSTVWY", freqs);
    set_ambiguity(log_space, 'X', "ACDEFGHIKLMNPQRSTVWY", freqs);
    break;

  case DNA_ALPH :
    set_ambiguity(log_space, 'R', "GA", freqs);
    set_ambiguity(log_space, 'Y', "TC", freqs);
    set_ambiguity(log_space, 'K', "GT", freqs);
    set_ambiguity(log_space, 'M', "AC", freqs);
    set_ambiguity(log_space, 'S', "GC", freqs);
    set_ambiguity(log_space, 'W', "AT", freqs);
    set_ambiguity(log_space, 'B', "GTC", freqs);
    set_ambiguity(log_space, 'D', "GAT", freqs);
    set_ambiguity(log_space, 'H', "ACT", freqs);
    set_ambiguity(log_space, 'V', "GCA", freqs);
    set_ambiguity(log_space, 'N', "ACGT", freqs);
    break;

  default :
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  }
  
  if (AMBIG_DEBUG) {
    fprintf(stderr, "Emission distribution: ");
    print_array(freqs, 5, 3, TRUE, stderr);
  }
}      

/********************************************************************
 * Take the counts from an ambiguous character and evenly distribute
 * them among the corresponding concreate characters.
 *
 * This function operates in log space.
 ********************************************************************/
void distribute_one_count
  (const char     ambig_char,
   const char*    concrete_chars,
   ARRAY_T* freqs)
{
  int    ambig_index;      /* Index of the ambiguous character. */
  PROB_T ambig_count;      /* Count of the ambiguous character. */
  int    num_concretes;    /* Number of associated definite characters. */
  int    i_concrete;       /* Which definite character are we considering? */
  char   concrete_char;    /* Current definite character. */
  int    concrete_index;   /* Index of the definite character. */
  PROB_T concrete_count;   /* Count of the definite character. */

  /* Get the count to be distributed. */
  ambig_index = alphabet_index(ambig_char, alphabet);
  ambig_count = get_array_item(ambig_index, freqs);

  /* Divide it by the number of corresponding concrete characters. */
  num_concretes = strlen(concrete_chars);
  ambig_count -= my_log2((PROB_T)num_concretes);

  /* Distribute it in equal portions to the given concrete characters. */
  for (i_concrete = 0; i_concrete < num_concretes; i_concrete++) {
    concrete_char = concrete_chars[i_concrete];
    concrete_index = alphabet_index(concrete_char, alphabet);
    concrete_count = get_array_item(concrete_index, freqs);

    /* Add the ambiguous counts. */
    concrete_count = LOG_SUM(concrete_count, ambig_count);

    set_array_item(concrete_index, concrete_count, freqs);
  }

  /* Set the ambiguous count to zero. */
  set_array_item(ambig_index, LOG_ZERO, freqs);
}
    


/********************************************************************
 * See .h file for description.
 ********************************************************************/
void distribute_ambiguous_counts
  (ARRAY_T* freqs)
{
  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    distribute_one_count('B', "DN", freqs);
    distribute_one_count('Z', "EQ", freqs);
    distribute_one_count('U', "ACDEFGHIKLMNPQRSTVWY", freqs);
    distribute_one_count('X', "ACDEFGHIKLMNPQRSTVWY", freqs);
    break;

  case DNA_ALPH :
    distribute_one_count('R', "GA", freqs);
    distribute_one_count('Y', "TC", freqs);
    distribute_one_count('K', "GT", freqs);
    distribute_one_count('M', "AC", freqs);
    distribute_one_count('S', "GC", freqs);
    distribute_one_count('W', "AT", freqs);
    distribute_one_count('B', "GTC", freqs);
    distribute_one_count('D', "GAT", freqs);
    distribute_one_count('H', "ACT", freqs);
    distribute_one_count('V', "GCA", freqs);
    distribute_one_count('N', "ACGT", freqs);
    break;

  default :
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  }
}
    


/********************************************************************
 * Set all the ambiguous characters to zero.
 ********************************************************************/
void zero_ambigs
  (ARRAY_T* freqs)
{
  int i_ambig;
  int num_chars = 0;
  int num_ambigs = 0;

  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    num_chars = NUM_AMINOS;
    num_ambigs = NUM_AMINO_AMBIGS;
    break;

  case DNA_ALPH :
    num_chars = NUM_BASES;
    num_ambigs = NUM_BASE_AMBIGS;
    break;

  default :
    carp(CARP_FATAL, "Alphabet uninitialized.\n");
  }
  
  for (i_ambig = 0; i_ambig < num_ambigs; i_ambig++) {
    set_array_item(num_chars + i_ambig, 0.0, freqs);
  }
}      

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
