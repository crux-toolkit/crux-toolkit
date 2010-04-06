/**
 * \file crux-utils.h
 * $Revision: 1.41 $
 * $Author: cegrant $
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
#include <algorithm>
#include "carp.h"
#include "utils.h"
#include "objects.h"
#include "peak.h"

#ifdef __cplusplus
  extern "C" {
#endif


/**
 *\returns a heap copy of the given string
 */
char* my_copy_string(const char* src);

/**
 * returns copy of the src string upto the specified length
 * includes a null terminating `\0' character
 * the string is heap allocated thus, user must free
 */
char* copy_string_part(const char* src, int length);

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * compare the absolute value of the difference of two numbers with an 
 * appropriate epsilon to get relations.
 * Multiplying the epsilon by sum of the comparands adjusts the comparison 
 * to the range of the numbers, allowing a single epsilon to be used for many, 
 * or perhaps all compares.
 */
int compare_float(FLOAT_T float_a, FLOAT_T float_b);

/**
 * Compares two numbers and returns TRUE if they are within the given
 * precision of each other, otherwise returns FALSE.  E.g. if
 * precision is 2, a and b must be equal when rounded to two decimal
 * places.
 */
BOOLEAN_T is_equal(FLOAT_T a, FLOAT_T b, int precision);

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
char** parse_filename_path(const char* file);

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
char** parse_filename_path_extension(const char* file, const char* extension);

/**
 * parses the filename
 * ex) ../../file_name => returns filename
 *\returns A heap allocated array of filename
 */
char* parse_filename(const char* file);

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
//void print_peptide_type(PEPTIDE_TYPE_T peptide_type, FILE* file);

/**
 * given two strings return a concatenated third string
 * \returns a heap allocated string that concatenates the two inputs
 */
char* cat_string(const char* string_one, const char* string_two);

/**
 * Adds the fileroot parameter to a string as a prefix.
 * Given a pointr to pointer to a string, if fileroot parameter is set
 * the memory for the string is reallocated, and the fileroot string
 * is added as a suffix.
 */
void prefix_fileroot_to_name(char** name);

/**
 * given the path and the filename return a file with path
 * "path/filename"
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(const char* path, const char* filename);

/**
 *\returns TRUE if float_a is between the interaval of min and max, else FALSE
 */
inline BOOLEAN_T compare_float_three(FLOAT_T float_a, FLOAT_T min, FLOAT_T max);

/**
 * returns the file size of the given filename
 */
long get_filesize(char *FileName);

/**
 * \brief A function for creating a directory to hold output files from crux.
 * 
 * Tries to create a the named directory for use as the output directory for crux.
 * If the overwrite option is true, an existing directory wtih that
 * name will not cause an error. 
 * 
 * \returns 0 if successful, -1 if an error occured.
*/
int create_output_directory(
  const char *output_folder, // Name of output folder.
  BOOLEAN_T overwrite  // Whether or not to overwrite an existing dir 
); 

/**
 * returns whether the given filename is a directory.
 * Returns TRUE if a directory, FALSE otherwise.
 * Terminates program if unable to determine status of file.
 */
BOOLEAN_T is_directory(const char *FileName);

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
  const char* fasta_filename,
  const char* name_tag,
  const char* file_extension,
  const char* suffix
  );

/**
 * \brief Take a filename, strip its leading path information (if
 * any) and file extension (if any).  Add a new path (if given) and a
 * new suffix (exension).
 *
 * If given ../dir/filename.ext, .new-ext, .ext, otherdir would return
 * otherdir/filename.new-ext 
 * \returns A heap allocated filename
 */
char* generate_name_path(
  const char* filename,
  const char* old_suffix,
  const char* new_suffix,
  const char* new_path
  );

/**
 * \brief Create the correct filename for a binary psm file, 
 * search.target.csm for target search and search.decoy-#.csm for 
 * decoy searches.
 *
 * Adds the appropriate
 * extension depending on the file index (0=target, 1=first decoy,
 * 2=second decoy, etc).
 * \returns A heap allocated char* with the new filename.
 */
char* generate_psm_filename(int file_index);

/**
 * \brief Open and create a file of the given name in the given
 * directory.
 *
 * Assumes the directory exists.  Fails if file can't be opened or if
 * file exists and overwrite is false.
 *\returns A file handle to the newly created file.
 */
FILE* create_file_in_path(
  const char* filename,  ///< the filename to create & open -in
  const char* directory,  ///< the directory to open the file in -in
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
  const char* string, ///< The string to compare -in
  const char* suffix  ///< The suffix to find in the string -in
  );

/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
BOOLEAN_T valid_peptide_sequence(const char* sequence);

/**
 * quickSort for FLOAT_Ts
 */
void quicksort(FLOAT_T numbers[], int array_size);

/**
 * \brief Shuffle an array of FLOAT_Ts.  Uses the Knuth algorithm.  Uses
 * get_random_number_interval() to generate random numbers. 
 */
void shuffle_floats(FLOAT_T* array, int size);

/**
 * \brief Comparison function for reverse sorting FLOAT_Ts.
 * \returns -1,0,1 if a is <,=,> b
 */
int compare_floats_descending(const void* a, const void* b);

/**
 *\returns a heap allocated feature name array for the algorithm type
 */
char** generate_feature_name_array();

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
    FLOAT_T* data, ///< the data to be fit. should be in descending order -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    FLOAT_T min_shift, ///< the minimum shift to allow -in
    FLOAT_T max_shift, ///< the maximum shift to allow -in
    FLOAT_T step, ///< the step for shift modification -in
    FLOAT_T corr_threshold, ///< minimum correlation, else no fit -in
    FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
    FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
    FLOAT_T* shift,     ///< the best shift -out
    FLOAT_T* correlation   ///< the best correlation -out
    );

/**
 * Fits a two-parameter Weibull distribution to the input data. 
 * \returns eta, beta and the correlation coefficient
 */
void fit_two_parameter_weibull(
    FLOAT_T* data, ///< the data to be fit -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    FLOAT_T shift, ///< the amount by which to shift our data -in
    FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
    FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
    FLOAT_T* correlation ///< the best correlation -out
    );

BOOLEAN_T string_to_mass_type(char*, MASS_TYPE_T*);
BOOLEAN_T mass_type_to_string(MASS_TYPE_T, char*);
//BOOLEAN_T string_to_peptide_type(char*, PEPTIDE_TYPE_T*);
//BOOLEAN_T peptide_type_to_string(PEPTIDE_TYPE_T type, char* type_str);
BOOLEAN_T string_to_sort_type(char*, SORT_TYPE_T*);
BOOLEAN_T sort_type_to_string(SORT_TYPE_T, char*);
BOOLEAN_T string_to_algorithm_type(char*, ALGORITHM_TYPE_T*);
BOOLEAN_T algorithm_type_to_string(ALGORITHM_TYPE_T, char*);
BOOLEAN_T string_to_scorer_type(char*, SCORER_TYPE_T*);
BOOLEAN_T scorer_type_to_string(SCORER_TYPE_T, char*);
BOOLEAN_T string_to_ion_type(char* , ION_TYPE_T*);
BOOLEAN_T ion_type_to_string(ION_TYPE_T, char*);

// new style of type_to_string and string_to_type functions
// requires an invalid value for each enum
char* command_type_to_file_string(COMMAND_T type);
const char* command_type_to_file_string_ptr(COMMAND_T type);
char* command_type_to_command_line_string(COMMAND_T type);
const char* command_type_to_command_line_string_ptr(COMMAND_T type);
COMMAND_T string_to_command_type(char*);
DIGEST_T string_to_digest_type(char*);
char* digest_type_to_string(DIGEST_T);
ENZYME_T string_to_enzyme_type(char*);
char* enzyme_type_to_string(ENZYME_T);
WINDOW_TYPE_T string_to_window_type(char*);
char* window_type_to_string(WINDOW_TYPE_T);


/**
 * \brief Open either the index or fasta file and prepare it for
 * searching.  Die if the input file cannot be found or read.
 * \returns The number of proteins in the file or index
 */
int prepare_protein_input(
  char* input_file,      ///< name of the fasta file or index directory
  INDEX_T** index,       ///< return new index here OR
  DATABASE_T** database);///< return new fasta database here

/**
 * \brief Perform the set-up steps common to all crux commands:
 * initialize parameters, parse command line, set verbosity, open
 * output directory, write params file. 
 *
 * Uses the given command name, arguments and options for parsing the
 * command line.
 */
void initialize_run(
  COMMAND_T cmd,              ///< the command we are initializing 
  const char** argument_list, ///< list of required arguments
  int num_arguments,          ///< number of elements in arguments_list
  const char** option_list,   ///< list of optional flags
  int num_options,            ///< number of elements in options_list
  int argc,                   ///< number of tokens on cmd line
  char** argv                 ///< array of command line tokens
);

/**
 *  Read the string of the form <first>-<last> and returns <first>
 *  or -1 if the range is invalid.
 */
int get_first_in_range_string(const char* range_string);
/**
 *  Read the string of the form <first>-<last> and returns <last>
 *  or -1 if the range is invalid.
 */
int get_last_in_range_string(const char* range_string);

/**
 * \brief  Decide if a spectrum has precursor charge of +1 or more (+2
 * or +3 or +4 etc). 
 * \returns 1 if spectrum precursor is singly charged or 0 if multiply
 * charged.
 */
int choose_charge(FLOAT_T precursor_mz, ///< m/z of spectrum precursor ion
                  PEAK_T* peaks,        ///< array of spectrum peaks
                  int num_peaks);       ///< size of peaks array

#ifdef __cplusplus
}
#endif

#endif
