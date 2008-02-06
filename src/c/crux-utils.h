/**
 * \file crux-utils.h
 * $Revision: 1.29 $
 * $Author: frewen $
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
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include "carp.h"
#include "utils.h"
#include "objects.h"

/**
 *\returns a heap copy of the given string
 */
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
 * \brief Parses the filename and path of given string.
 *  
 * The array returned, A, contains the filename (A[0]) and the path
 * (A[1]).  Path is NULL if no / in given name.  e.g. Given
 * "../../filname" returns A[0]="filename" and A[1]="../../".  Given
 * "filename" returns A[0] = "filename" and A[1] = NULL.
 *
 *\returns A heap allocated array of both filename and path
 */
char** parse_filename_path(char* file);

/**
 * \brief Parses the filename, path, and file extension of given string.
 *  
 * The array returned, A, contains the filename (A[0]) striped of the
 * given file extension and the path (A[1]).  Path is NULL if no / in
 * given name.  Filename is unchanged if it does not include the
 * extension. e.g. Given "../../filname.ext" and ".ext" returns
 * A[0]="filename" A[1]="../../".  Given "filename" returns A[0] =
 * NULL and A[1] = "filename". 
 *
 * \returns A heap allocated array of filename striped of extension and
 * path.
 */
char** parse_filename_path_extension(char* file, char* extension);

/**
 * parses the filename
 * ex) ../../file_name => returns filename
 *\returns A heap allocated array of filename
 */
char* parse_filename(char* file);

/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* int_to_char(unsigned int i);

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
  char* name_tag,
  char* file_extension,
  char* suffix
  );

/**
 * Open and create a file handle of a file that is named 
 * and located in user specified location
 * Assumes the directory exists
 *\returns a file handle of a file that is named and located in user specified location
 */
FILE* create_file_in_path(
  char* filename,  ///< the filename to create & open -in
  char* directory,  ///< the directory to open the file in -in
  BOOLEAN_T overwrite  ///< replace the file (T) or die if exists (F)
  );

/**
 * check if the string has the correct suffix
 * \returns TRUE, if the string starts with the suffix, else FALSE
 */
BOOLEAN_T prefix_compare(
  char* string, ///< The string to compare -in
  char* prefix  ///< The prefix to find in the string -in
  );

/**
 * check if the string has the correct suffix
 * \returns TRUE, if the string starts with the suffix, else FALSE
 */
BOOLEAN_T suffix_compare(
  char* string, ///< The string to compare -in
  char* suffix  ///< The suffix to find in the string -in
  );

/**
 * \brief Decide if a file name is a decoy csm file
 * \returns TRUE if name ends in -decoy-#.csm, else false
 */
BOOLEAN_T name_is_decoy(char* name);

/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
BOOLEAN_T valid_peptide_sequence( char* sequence);

/**
 * quickSort for floats
 */
void quicksort(float numbers[], int array_size);

/**
 *\returns a heap allocated feature name array for the algorithm type
 */
char** generate_feature_name_array(
  ALGORITHM_TYPE_T algorithm ///< the algorithm's feature name to produce -in
  );

/**
 *\returns the number of digits in the number
 */
int get_number_digits(
  int number ///< the number to count digits
  );

/**
 * User define our upper and our lower bounds.
 * The random number will always be 
 * between low and high, inclusive.
 * There is no seeding in this function, user must do it for themselves
 *\returns a random number between the interval user provides
 */
int get_random_number_interval(
  int low, ///< the number for lower bound -in
  int high ///< the number for higher bound -in
  );

/**
 * Fits a three-parameter Weibull distribution to the input data. 
 * \returns eta, beta, c (which in this case is the amount the data should
 * be shifted by) and the best correlation coefficient
 */

void fit_three_parameter_weibull(
    float* data, ///< the data to be fit. should be in descending order -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    float min_shift, ///< the minimum shift to allow -in
    float max_shift, ///< the maximum shift to allow -in
    float step, ///< the step for shift modification -in
    float* eta,      ///< the eta parameter of the Weibull dist -out
    float* beta,      ///< the beta parameter of the Weibull dist -out
    float* shift,     ///< the best shift -out
    float* correlation   ///< the best correlation -out
    );

/**
 * Fits a two-parameter Weibull distribution to the input data. 
 * \returns eta, beta and the correlation coefficient
 */
void fit_two_parameter_weibull(
    float* data, ///< the data to be fit -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    float shift, ///< the amount by which to shift our data -in
    float* eta,      ///< the eta parameter of the Weibull dist -out
    float* beta,      ///< the beta parameter of the Weibull dist -out
    float* correlation ///< the best correlation -out
    );

BOOLEAN_T string_to_mass_type(char*, MASS_TYPE_T*);
BOOLEAN_T mass_type_to_string(MASS_TYPE_T, char*);
BOOLEAN_T string_to_peptide_type(char*, PEPTIDE_TYPE_T*);
BOOLEAN_T peptide_type_to_string(PEPTIDE_TYPE_T type, char* type_str);
BOOLEAN_T string_to_sort_type(char*, SORT_TYPE_T*);
BOOLEAN_T sort_type_to_string(SORT_TYPE_T, char*);
BOOLEAN_T string_to_algorithm_type(char*, ALGORITHM_TYPE_T*);
BOOLEAN_T algorithm_type_to_string(ALGORITHM_TYPE_T, char*);
BOOLEAN_T string_to_scorer_type(char*, SCORER_TYPE_T*);
BOOLEAN_T scorer_type_to_string(SCORER_TYPE_T, char*);
BOOLEAN_T string_to_output_type(char*, MATCH_SEARCH_OUTPUT_MODE_T*);
BOOLEAN_T output_type_to_string(MATCH_SEARCH_OUTPUT_MODE_T, char*);
BOOLEAN_T string_to_ion_type(char* , ION_TYPE_T*);
BOOLEAN_T ion_type_to_string(ION_TYPE_T, char*);

#endif
