/**
 * \file crux-utils.h
 * $Revision: 1.15 $
 * $Author: cpark $
 * \brief Utilities for the crux project
 */
#ifndef CRUX_UTILS_H
#define CRUX_UTILS_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "objects.h"

char* my_copy_string(char* src);


/**
 * returns copy of the src string upto the specified length
 * includes a null terminating `\0' character
 * the string is heap allocated thus, user must free
 */
char* copy_string_part(char*src, int length);

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * compare the absolute value of the difference of two numbers with an 
 * appropriate epsilon to get relations.
 * Multiplying the epsilon by sum of the comparands adjusts the comparison 
 * to the range of the numbers, allowing a single epsilon to be used for many, 
 * or perhaps all compares.
 */
int compare_float(float float_a, float float_b);

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * fast and simple, but some limitations. Assumes,
 * "Two floats in memory, interpret their bit pattern as integers, 
 * and compare them, we can tell which is larger"
 */
inline int compare_float_fast(float float_a, float float_b);

/**
 * parses the filename and path  
 * returns an array A, with A[0] the filename and A[1] the path to the filename
 * returns A[1] NULL if only a filename was passed in
 * ex) ../../file_name => returns filename , ../../
 *     file_name => returns filename, NULL
 *\returns A heap allocated array of both filename and path
 */
char** parse_filename_path(char* file);

/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* int_to_char(int i);

/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* signed_int_to_char(int i);

/**
 *prints the peptide type given it's enum value
 */
void print_peptide_type(PEPTIDE_TYPE_T peptide_type, FILE* file);

/**
 * given two strings return a concatenated third string
 * \returns a heap allocated string that concatenates the two inputs
 */
char* cat_string(char* string_one, char* string_two);

/**
 * given the path and the filename return a file with path
 * "path/filename"
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(char* path, char* filename);

/**
 *\returns TRUE if float_a is between the interaval of min and max, else FALSE
 */
inline BOOLEAN_T compare_float_three(float float_a, float min, float max);

/**
 * returns the file size of the given filename
 */
long get_filesize(char *FileName);

/**
 * deletes a given directory and it's files inside.
 * assumes that there's no sub directories, only files
 * \returns TRUE if successfully deleted directory
 */
BOOLEAN_T delete_dir(char* dir);

/**
 * given a fasta_file name it returns a name with the name_tag add to the end
 * format: myfasta_nameTag
 * \returns A heap allocated file name of the given fasta file
 */
char* generate_name(
  char* fasta_filename,
  char* name_tag
  );

/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
BOOLEAN_T valid_peptide_sequence( char* sequence);

/**
 *
 *quickSort for floats
 */
void quicksort(float numbers[], int array_size);
#endif
