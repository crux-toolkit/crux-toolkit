/**
 * \file utils.h
 * \brief Various useful generic utilities.
 ********************************************************************/
// AUTHOR: William Stafford Noble
// CREATE DATE: 9-8-97

#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#ifndef linux
#include <ieeefp.h>
#endif*/

// Macro allowing us to build using floats or double
#ifdef USE_DOUBLES
typedef double FLOAT_T;
#else
typedef float FLOAT_T;
#endif

#define FALSE 0
#define TRUE 1
typedef short BOOLEAN_T;

typedef int VERBOSE_T;
#define INVALID_VERBOSE 0
#define QUIET_VERBOSE 1
#define NORMAL_VERBOSE 2
#define HIGH_VERBOSE 3
#define HIGHER_VERBOSE 4
#define DUMP_VERBOSE 5

#define FILENAME_LENGTH 4096
#define BAD_SCORE -1
#define IDLENGTH 256
#define PEPTIDELENGTH 80
#define LINELENGTH 4096

#define EPS 1E-10

#define BIG HUGE_VAL
#define LITTLE -BIG

extern int verbosity;

#ifdef DARWIN
/*********************************************************
 This function replaces the GNU extension of the same name.
 Reads a line from the given stream.
 *********************************************************/
int getline(char **lineptr, size_t *n, FILE *stream);
#endif

/***********************************************************************
 * Return a not-a-number.
 ***********************************************************************/
double NaN
  (void);

/***********************************************************************
 * Return elapsed time in microseconds since the last call.
 ***********************************************************************/
double wall_clock(void);

/************************************************************************
 * int open_file
 *
 * Open a file gracefully.
 *
 * RETURN: Was the open successful?
 ************************************************************************/
BOOLEAN_T open_file
(char*     filename,            // Name of the file to be opened.
 char*     file_mode,           // Mode to be passed to fopen.
 BOOLEAN_T allow_stdin,         // If true, filename "-" is stdin.
 char*     file_description,   
 char*     content_description,
 FILE**    afile);              // Pointer to the open file.

/********************************************************************
 * DEBUG_CODE (macro)
 *
 * Allow debugging code to be included or excluded from a compiled
 * program.
 ********************************************************************/
#ifdef DEBUG
#define DEBUG_CODE( debug_value, code_fragment ) \
   { if (debug_value) { code_fragment } }
#else
#define DEBUG_CODE( debug_value, code_fragment )
#endif 

/********************************************************************
 * Allocate dynamic memory. Die gracefully if memory is exhausted.
 ********************************************************************/
void *mymalloc
  (size_t size);
void *mycalloc
  (size_t nelem,
   size_t size);
void * myrealloc
  (void * ptr,
   size_t size);

/********************************************************************
 * fwrite with a check to make sure it was successful (useful for NFS problems)
 ********************************************************************/
BOOLEAN_T myfwrite
  (const void *ptr, 
   size_t size, 
   size_t nitems, 
   FILE *stream);

/***************************************************************************
 * Dynamically create or grow an array; 
 * P = pointer, N = new size, T = type 
 **************************************************************************/
typedef void *malloc_t;
#define Resize(P,N,T) { \
  void *new_P; \
  new_P = (P) ? realloc((malloc_t)(P), (N)*sizeof(T)) : malloc((N)*sizeof(T)); \
  if (!new_P) { \
    fprintf(stderr, "Resize(" #P "," #N "," #T ") failed!\n"); \
    exit(1); \
  } \
  (P) = (T *) new_P; \
}

/********************************************************************
 * Only free memory if the given pointer is non-null.
 ********************************************************************/
#define myfree(x) if (x) free((char* ) (x))

/********************************************************************
 * Set the seed for the random number generator.
 ********************************************************************/
void my_srand
  (long seed);

/********************************************************************
 * Get a random number X such that 0 <= X < 1.
 ********************************************************************/
double my_drand
  (void);

/********************************************************************
 * Math macros.
 ********************************************************************/
/* Note that the following type must be the  same as the MTYPE and ATYPE
   defined in 'matrix.h' and 'array.h'. */
typedef double PROB_T;       // Type definition for probability/frequency.
#define PROB_SCAN " %lf"     // Scanf string for PROB_T.

#define LOG_ZERO  (-1.0E10)  // Zero on the log scale.
#define LOG_SMALL (-0.5E10)  // Threshold below which everything is zero.
#define BITS      (33.2)     // = LOG2(-LOG_ZERO)

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/***************************************************************************
 * Find the nearest integer.
 ***************************************************************************/
#define nint(x) ((int)((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)))

/**************************************************************************
 * Compute the logarithm of x, when 0 <= x <= 1.
 **************************************************************************/
void init_log_prob
  (void);

PROB_T log_prob
  (PROB_T value);

#define log_prob2(x)   (log_prob(x) * 1.44269504)

/**************************************************************************
 * Compute the logarithm of x.  Returns LOG_ZERO if x==0.
 **************************************************************************/
PROB_T my_log
  (PROB_T x);

#define my_log2(x)   (my_log(x) * 1.44269504)

#define EXP2(x)                          \
( ( (x) < LOG_SMALL) ?                   \
  0.0 :                                  \
  (exp((x) * 0.69314718 ))               \
) 

/**************************************************************************
 * Given the logs (in base 2) of two numbers, return the log of their
 * sum.
 *
 * This function is optimized based upon the following formula:
 *
 *      log(x+y) = log(x) + log(1 + exp(log(y) - log(x)))
 *
 **************************************************************************/
#define LOG_VALUE(logx) \
( ( (logx) < LOG_SMALL ) ? \
    LOG_ZERO : \
    (logx) \
)

#define LOG_SUM1(logx, logy) \
( \
  ( ( (logx) - (logy) ) > BITS ) ? \
    LOG_VALUE(logx) : \
    (logx) + my_log2( 1 + EXP2((logy) - (logx) ) ) \
)

#define LOG_SUM(logx, logy) \
( \
  ( (logx) > (logy) ) ? \
    LOG_SUM1( (logx), (logy) ) : \
    LOG_SUM1( (logy), (logx) ) \
)

/**************************************************************************
 * Test for zero on a value that may be either a log or a raw float.
 **************************************************************************/
BOOLEAN_T is_zero
  (double    value,
   BOOLEAN_T log_form);

/**************************************************************************
 * Test to see if two values are approximately equal.
 **************************************************************************/
BOOLEAN_T almost_equal
  (double value1,
   double value2,
   double slop);

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
char*  boolean_to_string
 (BOOLEAN_T the_boolean);

BOOLEAN_T boolean_from_string
  (char* true_or_false);

/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
BOOLEAN_T char_in_string
  (const char* a_string,
   char        a_char);

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
const char*  convert_enum_type
  (int     enum_type,  /* The enumerated type object to be converted. */
   const char*  enum_strs[],  /* String values associated with this type. */
   int     num_enums); /* Number of values of the type. */

int convert_enum_type_str
  (const char*   enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   const char**  enum_strs,     /* String values associated with this type. */
   int     num_enums);    /* Number of values of the type. */

/**************************************************************************
 * Get the name of the CPU.
 **************************************************************************/
const char* hostname(void);

/**************************************************************************
 * Get the current date and time.
 **************************************************************************/
const char* date_and_time(void);

/**************************************************************************
 * Copy a string, with allocation.
 **************************************************************************/
void copy_string
 (char**  target,
  char*   source);

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target);

/************************************************************************
 * parses a file of length max_lines and returns an array of strings
 ************************************************************************/
char** parse_file(
  char* file_name,
  int max_lines, 
  int* num_lines
  );

#ifdef __cplusplus
}
#endif

#endif

