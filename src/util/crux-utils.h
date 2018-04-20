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
#ifndef _MSC_VER
#include <dirent.h>
#include <unistd.h>
#endif
#include <time.h>
#include <algorithm>
#include <limits>
#include "io/carp.h"
#include "utils.h"
#include "model/objects.h"
#include "parameter.h"
#include "StringUtils.h"
#include "model/Peak.h"
#include "app/CruxApplication.h"

#include <set>
#include <sstream>
#include <string>
#include <vector>

/**
 * The number of features used to represent a PSM for Percolator or q-ranker.
 */
const unsigned int NUM_FEATURES = 16;

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
bool is_equal(FLOAT_T a, FLOAT_T b, int precision);

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
char** parse_filename_path(const std::string& file);

/**
 * \brief Parses the filename, path, and file extension of given string.
 *  
 * The array returned, A, contains the filename (A[0]) striped of the
 * given file extension and the path (A[1]).  If extension is NULL,
 * strips all characters after the last "." from the filename.  Use
 * parse_filename_path() to return filename with extension.  
 * Returned path is NULL if no "/" in given name.
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
 * Examines filename to see if it ends in the given extension
 * \returns True if filename ends in the extension, else false.
 */
bool has_extension(std::string filename, std::string extension);

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
 */
std::string prefix_fileroot_to_name(const std::string& name);

/**
 * \returns the filepath 'output_dir'/'fileroot'.'filename' 
 */
std::string make_file_path(
  const std::string& filename ///< the name of the file
  );


/**
 * given the path and the filename return a file with path
 * "path/filename"
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(const char* path, const char* filename);

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
  const std::string& output_folder, // Name of output folder.
  bool overwrite  // Whether or not to overwrite an existing dir 
); 

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
 * any) and file extension (if any).  Tries all file extensions until
 * one is found.  Add a new path (if given) and a new suffix (exension).
 *
 * If given ../dir/filename.ext, [.txt, .ext, t], .new-ext, otherdir
 * would return  otherdir/filename.new-ext 
 * \returns A heap allocated filename
 */
char* generate_name_path(
  const char* filename,
  std::vector<const char*> old_suffixes,
  const char* new_suffix,
  const char* new_path
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
 * \brief Open and create a file of the given name in the given
 * directory.
 *
 * Assumes the directory exists.  Fails if file can't be opened or if
 * file exists and overwrite is false.
 *\returns A file handle to the newly created file.
 */
FILE* create_file_in_path(
  const std::string& filename,  ///< the filename to create & open -in
  const std::string& directory,  ///< the directory to open the file in -in
  bool overwrite  ///< replace the file (T) or die if exists (F)
  );

/**
 * \brief c++ version of create_file_in_path
 */
std::ofstream* create_stream_in_path(
  const char* filename,  ///< the filename to create & open -in
  const char* directory,  ///< the directory to open the file in -in
  bool overwrite  ///< replace file (T) or die if exists (F)
  );

/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
bool valid_peptide_sequence(const std::string& sequence);

/**
 * quickSort for FLOAT_Ts
 */
void quicksort(FLOAT_T numbers[], int array_size);

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
 * \brief Shuffle an array of FLOAT_Ts.  Uses the Knuth algorithm.  Uses
 * get_random_number_interval() to generate random numbers. 
 */
void shuffle_floats(FLOAT_T* array, int size);

/**
 * \brief Shuffles an array of elements.  Uses the Knuth algorithm.  Uses
 * get_random_number_interval() to generate random numbers. 
 */
template<typename T>
void shuffle_array(T* array, int size) {
  if (array == NULL) {
    carp(CARP_ERROR, "Cannot shuffle NULL array.");
    return;
  }

  int idx, switch_idx;
  int last_element_idx = size - 1;
  T temp_value;
  for (idx = 0; idx < size; idx++) {
    switch_idx = get_random_number_interval(idx, last_element_idx);
    temp_value = array[idx];
    array[idx] = array[switch_idx];
    array[switch_idx] = temp_value;
  }
}

/**
 *\returns the number of digits in the number
 */
int get_number_digits(
  int number ///< the number to count digits
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

bool string_to_mass_type(const std::string& name, MASS_TYPE_T*);
bool mass_type_to_string(MASS_TYPE_T, char*);
bool string_to_algorithm_type(char*, ALGORITHM_TYPE_T*);
bool algorithm_type_to_string(ALGORITHM_TYPE_T, char*);
bool string_to_scorer_type(const std::string& name, SCORER_TYPE_T*);
const char* scorer_type_to_string(SCORER_TYPE_T);
bool string_to_ion_type(const std::string& name, ION_TYPE_T*);
bool ion_type_to_string(ION_TYPE_T, char*);
char* ion_type_to_string(ION_TYPE_T type);

// new style of type_to_string and string_to_type functions
// requires an invalid value for each enum
DIGEST_T string_to_digest_type(const std::string& name);
const char* digest_type_to_string(DIGEST_T);
ENZYME_T string_to_enzyme_type(const std::string& name);
const char* enzyme_type_to_string(ENZYME_T);
OBSERVED_PREPROCESS_STEP_T string_to_observed_preprocess_step(const std::string& name);
char* observed_preprocess_step_to_string(OBSERVED_PREPROCESS_STEP_T type);
WINDOW_TYPE_T string_to_window_type(const std::string& name);
PARSIMONY_TYPE_T string_to_parsimony_type(const std::string& name);
MEASURE_TYPE_T string_to_measure_type(const std::string& name);
char * measure_type_to_string(MEASURE_TYPE_T type);
THRESHOLD_T string_to_threshold_type(const std::string& name);
char * threshold_type_to_string(THRESHOLD_T type);
QUANT_LEVEL_TYPE_T string_to_quant_level_type(const std::string& name);
COLTYPE_T string_to_column_type(const std::string& name);
COMPARISON_T string_to_comparison(const std::string& name);
DECOY_TYPE_T string_to_decoy_type(const std::string& name);
DECOY_TYPE_T string_to_tide_decoy_type(const std::string& name);
char* decoy_type_to_string(DECOY_TYPE_T type);
MASS_FORMAT_T string_to_mass_format(const std::string& name);
char* mass_format_type_to_string(MASS_FORMAT_T type);
SCORE_FUNCTION_T string_to_score_function_type(const std::string& name);
char* score_function_type_to_string(SCORE_FUNCTION_T type);

HARDKLOR_ALGORITHM_T string_to_hardklor_algorithm_type(const std::string& name);
std::string hardklor_hardklor_algorithm_type_to_string(HARDKLOR_ALGORITHM_T type);

/**
 * \brief Open the fasta file and prepare it for
 * searching.  Die if the input file cannot be found or read.
 * \returns The number of proteins in the file
 */
int prepare_protein_input(
  const std::string& input_file,      ///< name of the fasta file
  Database** database);///< return new fasta database here

/**
 * converts a string in #-# format to
 * a first and last variable
 * \returns whether the extraction was successful or not
 */

bool get_first_last_scan_from_string(
  const std::string& const_scans_string,
  int& first_scan,
  int& last_scan
  );

bool get_scans_from_string(
  const std::string& const_scans_string,
  std::set<int>& scans
);


template<typename TValue>
static bool get_range_from_string(
  const std::string& const_range_string, ///< the string to extract 
  TValue& first,  ///< the first value
  TValue& last ///< the last value
  ) {

  if (const_range_string.empty()) {
    first = (TValue)0;
    last = std::numeric_limits<TValue>::max();
    return true;
  }
  char* range_string = my_copy_string(const_range_string.c_str());

  bool ret;

  char* dash = strchr(range_string, '-');
  if (dash == NULL) { // a single number
    ret = StringUtils::TryFromString(range_string, &first);
    last = first;
  } else {
    //invalid if more than one dash
    const char* dash_check = strchr(dash + 1, '-');
    if (dash_check) {
      ret = false;
    } else {
      *dash = '\0';
      ret = StringUtils::TryFromString(range_string, &first);
      *dash = '-';
      dash++;
      ret &= StringUtils::TryFromString(dash, &last);
    }
  }

  free(range_string);
  return ret;    
}

/**
 *\brief Extend a given string with lines not exceeding a specified width, 
 * breaking on spaces.
 */
void strcat_formatted(
  char*       string_to_extend,
  const char* lead_string,        // Appears at the start of each line.
  const char* extension           // Text to add.
);

/**
 * \brief Checks if the given input file contains target, decoy PSMs or 
 * concatenated search results.
 *
 *\returns corrected file names. It does not check if files are exist.
 */
void check_target_decoy_files(
  std::string &target,   //filename of the target PSMs
  std::string &decoy     //filename of the decoy PSMs
);

void get_search_result_paths(
  const std::vector<std::string>& infiles,
  std::vector<std::string>& outpaths ///< paths of all search results -out
);


/**
 * \brief Checks if the given input file contains target, decoy PSMs or 
 * concatenated search results.
 *
 *\returns corrercted file names. It does not check if files are exist.
 */
void check_target_decoy_files(
  std::string &target,   //filename of the target PSMs
  std::string &decoy     //filename of the decoy PSMs
);

void get_files_from_list(
  const std::string &infile, ///< path of the first file.
  std::vector<std::string> &outpaths ///< paths of all search results -out
);

bool parseUrl(std::string url, std::string* host, std::string* path);
std::string httpRequest(const std::string& url, const std::string& data = "", bool waitForResponse = true);
void postToAnalytics(const std::string& appName);

#endif
