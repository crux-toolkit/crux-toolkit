/**
 * \file crux-utils.cpp
 * \brief General-use functions for crux
 */

#include <sys/types.h>
#include <fstream>
#include <errno.h>
#include <sys/stat.h>
#ifndef _MSC_VER
#include <unistd.h>
#include <sys/dir.h>
#include <dirent.h>
#endif
#include <iostream>
#include "crux-utils.h"
#include "model/Database.h"
#include "Params.h"
#include "parameter.h"
#include "Params.h"
#include "StringUtils.h"
#include "WinCrux.h"
#include "io/LineFileReader.h"
#define  NO_BOOST_DATE_TIME_INLINE
#include <boost/filesystem.hpp>
#include <boost/asio.hpp>
#include "FileUtils.h"
#include "crux_version.h"

using namespace std;

/**
 * PRECISION, determines the precision of the compare float, users
 * should lower the number if need more precision
 */
static const FLOAT_T PRECISION = 0.000000005; 

/**
 * the maximum error in terms of Units in the Last Place. 
 * This specifies how big an error we are willing to accept in terms of the value of the least significant 
 * digit of the floating point numbers representation. 
 * MAX_ULPS can also be interpreted in terms of how many representable floats 
 * we are willing to accept between A and B. This function will allow MAX_ULPS-1 floats between A and B.
 */
static const int MAX_ULPS = 2;

static const unsigned int TARGET_STRING_LENGTH = 6; //The length of string "target"
static const unsigned int DECOY_STRING_LENGTH = 5; // The length of string "decoy"                                                                                                                                                                                             


/* Functions for converting custom types to and from strings */

static const int INVALID_ENUM_STRING = -10;
/**
 * The string version of the score functions
 * Added by Andy Lin
 */

static const char* score_function_strings[NUMBER_SCORE_FUNCTIONS] = {
  "invalid", "xcorr", "residue-evidence", "both"
};

SCORE_FUNCTION_T string_to_score_function_type(const string& name) {
  int score_function_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
                                                 score_function_strings,
                                                 NUMBER_SCORE_FUNCTIONS);
  if (score_function_int < 0) {
    score_function_int = 0;
  }
  return (SCORE_FUNCTION_T)score_function_int;
}

char* score_function_type_to_string(SCORE_FUNCTION_T type) {
  return my_copy_string(score_function_strings[type]);
}

/**
 * The string version of the decoy types
 */
static const char* decoy_type_strings[NUMBER_DECOY_TYPES] = {
  "invalid", "none", "reverse", "protein-shuffle",
  "peptide-shuffle", "peptide-reverse"
};

DECOY_TYPE_T string_to_decoy_type(const string& name) {
  int decoy_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
                                        decoy_type_strings,
                                        NUMBER_DECOY_TYPES);
  if (decoy_int < 0) {
    decoy_int = 0;
  }

  return (DECOY_TYPE_T)decoy_int;
}

DECOY_TYPE_T string_to_tide_decoy_type(const string& name) {
  if (name == "none") {
    return NO_DECOYS;
  } else if (name == "shuffle") {
    return PEPTIDE_SHUFFLE_DECOYS;
  } else if (name == "peptide-reverse") {
    return PEPTIDE_REVERSE_DECOYS;
  } else if (name == "protein-reverse") {
    return PROTEIN_REVERSE_DECOYS;
  }
  carp(CARP_FATAL, "Invalid decoy type %s", name.c_str());
  return(INVALID_DECOY_TYPE); // Avoid compiler warning.
}

char* decoy_type_to_string(DECOY_TYPE_T type) {
  return my_copy_string(decoy_type_strings[type]);
}

/**
 * The string version of the mass format types
 */
static const char* mass_format_type_strings[NUMBER_MASS_FORMATS] = 
  { "invalid", "mod-only", "total", "separate" };

MASS_FORMAT_T string_to_mass_format(const string& name) {
  int mass_format_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                              mass_format_type_strings, 
                                              NUMBER_MASS_FORMATS);
  if (mass_format_int < 0) {
    mass_format_int = 0;
  }

  return (MASS_FORMAT_T)mass_format_int;
}

char* mass_format_type_to_string(MASS_FORMAT_T type) {
  return my_copy_string(mass_format_type_strings[type]);
}

/**
 * The string version of isotopic mass type (average, mono)
 */
static const char* mass_type_strings[NUMBER_MASS_TYPES] = {"average", "mono"};

bool string_to_mass_type(const string& name, MASS_TYPE_T* result) {
  bool success = true;
  //this is copied from parameter.c::get_peptide_mass_type
  int mass_type = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                        mass_type_strings, NUMBER_MASS_TYPES);

  (*result) = (MASS_TYPE_T)mass_type;

  if (mass_type < 0) {
    success = false;
  }
  return success;
}

bool mass_type_to_string(MASS_TYPE_T type, char* type_str) {
  bool success = true;
  if ((int)type > NUMBER_MASS_TYPES) {
    success = false;
    type_str = NULL;
  }

  //  type_str = mass_type_strings[type];
  strcpy(type_str, mass_type_strings[type]);

  return success;
}

// replace peptide type
/**
 * The string versions of digest types
 */
static const char* digest_type_strings[NUMBER_DIGEST_TYPES] =
  {"invalid", "full-digest", "partial-digest", "non-specific-digest"};

DIGEST_T string_to_digest_type(const string& name) {
  int clev_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                       digest_type_strings, 
                                       NUMBER_DIGEST_TYPES);
  if (clev_int < 0) {
    clev_int = 0;
  }

  return (DIGEST_T)clev_int;
}

const char* digest_type_to_string(DIGEST_T type){
  if( (int)type > NUMBER_DIGEST_TYPES){
    return NULL;
  }

  const  char* type_str = digest_type_strings[type];
  return type_str;
}

/**
 * The string version of observed preprocessing types
 */
static const char* observed_preprocess_step_strings[NUMBER_PREPROCESS_STEPS] =
  {"invalid", "discretize", "remove-precursor", "square-root", "remove-grass",
   "ten-bin", "xcorr"};

OBSERVED_PREPROCESS_STEP_T string_to_observed_preprocess_step(const string& name) {

  int obs_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
				      observed_preprocess_step_strings,
				      NUMBER_PREPROCESS_STEPS);
  if ( obs_int < 0 ){
    obs_int = 0;
  }
  return (OBSERVED_PREPROCESS_STEP_T)obs_int;
}

char* observed_preproces_step_to_string(OBSERVED_PREPROCESS_STEP_T type) {
  if ( (int)type > NUMBER_PREPROCESS_STEPS) {
    return NULL;
  }
  char* obs_str = my_copy_string(observed_preprocess_step_strings[type]);
  return(obs_str);
}

/**
 * The string version of enzyme types
 */
static const char* enzyme_type_strings[NUMBER_ENZYME_TYPES] = {
  "invalid", "no-enzyme", "trypsin", "trypsin/p", "chymotrypsin",
  "elastase", "clostripain", "cyanogen-bromide", "iodosobenzoate",
  "proline-endopeptidase", "staph-protease", "asp-n", "lys-c",
  "lys-n" , "arg-c" , "glu-c" , "pepsin-a", "elastase-trypsin-chymotrypsin",
  "custom-enzyme"
};

ENZYME_T string_to_enzyme_type(const string& name) {
  int enz_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                      enzyme_type_strings, 
                                      NUMBER_ENZYME_TYPES);
  if (enz_int < 0) {
    enz_int = 0;
  }

  return (ENZYME_T)enz_int;
}

const char* enzyme_type_to_string(ENZYME_T type){
  if ((int)type > NUMBER_ENZYME_TYPES) {
    return NULL;
  }
  const char* type_str = enzyme_type_strings[type];
  return type_str;
}

/**
 * The string version of mass types
 */
static const char* window_type_strings[NUMBER_WINDOW_TYPES] = 
  {"invalid", "mass", "mz", "ppm"};

WINDOW_TYPE_T string_to_window_type(const string& name) {
  int window_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                         window_type_strings, 
                                         NUMBER_WINDOW_TYPES);
  if (window_int < 0) {
    window_int = 0;
  }

  return (WINDOW_TYPE_T)window_int;
}

/**
 * The string version of parsimony types
 */
static const char* parsimony_type_strings[NUMBER_PARSIMONY_TYPES] =
  {"invalid", "simple", "greedy", "none"};

PARSIMONY_TYPE_T string_to_parsimony_type(const string& name) {
  int parsimony_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
                                            parsimony_type_strings,
                                            NUMBER_PARSIMONY_TYPES);
  if (parsimony_int < 0) {
    parsimony_int = 0;
  }

  return (PARSIMONY_TYPE_T)parsimony_int;
}

/**
 * The string version of measure types
 */
static const char* measure_type_strings[NUMBER_MEASURE_TYPES] =
  {"invalid", "RAW", "SIN", "NSAF", "dNSAF", "EMPAI"};

MEASURE_TYPE_T string_to_measure_type(const string& name) {
  int measure_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
                                          measure_type_strings,
                                          NUMBER_MEASURE_TYPES);
  if (measure_int < 0) {
    measure_int = 0;
  }
  
  return (MEASURE_TYPE_T)measure_int;
}


char * measure_type_to_string(MEASURE_TYPE_T type) {
  if ((int)type > NUMBER_MEASURE_TYPES) {
    return NULL;
  }

  char * type_str = my_copy_string(measure_type_strings[type]);
  
  return type_str;
}

/**
 * the string version of threshold type
 */
static const char* threshold_type_strings[NUMBER_THRESHOLD_TYPES] = {
  "invalid", "none", "qvalue", "custom"
};

THRESHOLD_T string_to_threshold_type(const string& name) {

  int threshold_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
    threshold_type_strings,
    NUMBER_THRESHOLD_TYPES);

  if ( threshold_int < 0) {
    threshold_int = 0;
  }
  return (THRESHOLD_T)threshold_int;

}

char* threshold_type_to_string(THRESHOLD_T type) {
  if ( (int)type > NUMBER_THRESHOLD_TYPES) {
    return NULL;
  }

  char* type_str = my_copy_string(threshold_type_strings[type]);

  return type_str;
}

/**
 * The string version of quantification level  types
 */
static const char* quant_level_type_strings[NUMBER_QUANT_LEVEL_TYPES] =
  {"invalid", "peptide", "protein"};

QUANT_LEVEL_TYPE_T string_to_quant_level_type(const string& name) {
  int quant_int = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
                                        quant_level_type_strings,
                                        NUMBER_QUANT_LEVEL_TYPES);
  if (quant_int < 0) {
    quant_int = 0;
  }
  return (QUANT_LEVEL_TYPE_T)quant_int;
}

/**
 * The string version of column types
 */ 
static const char* column_type_strings[NUMBER_COLTYPES] =
  {"invalid", "int", "real", "string"};

COLTYPE_T string_to_column_type(const string& name) {
  int coltype_int = convert_enum_type_str(
    name.c_str(),
    INVALID_ENUM_STRING,
    column_type_strings,
    NUMBER_COLTYPES);

  if (coltype_int < 0) {
    coltype_int = 0;
  }

  return (COLTYPE_T)coltype_int;

}

/**
 * The string version of comparison types
 */
static const char* comparison_type_strings[NUMBER_COMPARISONS] =
  {"invalid", "lt", "lte", "eq", "gte", "gt", "neq"};

COMPARISON_T string_to_comparison(const string& name) {
  int comparison_int = convert_enum_type_str(
    name.c_str(),
    INVALID_ENUM_STRING,
    comparison_type_strings,
    NUMBER_COMPARISONS);

  if (comparison_int < 0) {
    comparison_int = 0;
  }

  return (COMPARISON_T)comparison_int;
  
}

/*
 * The string version of ion types
 */
static const char* ion_type_strings[NUMBER_ION_TYPES] = {
  "a", "b", "c", "x", "y", "z", "p", "by", "bya", "all" };

bool string_to_ion_type(const string& name, ION_TYPE_T* result) {
  bool success = true;
  int ion_type = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                       ion_type_strings, NUMBER_ION_TYPES);
  (*result) = (ION_TYPE_T)ion_type;
  if (ion_type < 0) {
    success = false;
  }
  return success;
}

bool ion_type_to_string(ION_TYPE_T type,
                        char* type_str) {
  bool success = true;
  if ((int)type > NUMBER_ION_TYPES) {
    success = false;
    type_str = NULL;
  }
  strcpy(type_str, ion_type_strings[type]);
  return success;
}

/*
 * The string version of ALGORITHM_TYPE_T
 */
static const char* algorithm_type_strings[NUMBER_ALGORITHM_TYPES] = {
  "percolator", "rczar", "curve-fit",
  "none", "all", "q-ranker"
};

bool string_to_algorithm_type(char* name, ALGORITHM_TYPE_T* result) {
  bool success = true;

  int algorithm_type = convert_enum_type_str(name, INVALID_ENUM_STRING,
                                             algorithm_type_strings,
                                             NUMBER_ALGORITHM_TYPES);
  *result = (ALGORITHM_TYPE_T)algorithm_type;
  if (algorithm_type < 0) {
    success = false;
  }
  return success;
}

bool algorithm_type_to_string(ALGORITHM_TYPE_T type, char* type_str) {
  bool success = true;
  if ((int)type > NUMBER_ALGORITHM_TYPES) {
    success = false;
    type_str = NULL;
  }
  strcpy(type_str, algorithm_type_strings[type]);
  return success;
}

/* 
 * The string version of HARDKLOR_ALGORITHM_T
 */
static const char* hardklor_algorithm_type_strings[NUMBER_HK_ALGORITHM_TYPES] = {
  "invalid", "basic", "fewest-peptides", "fast-fewest-peptides", 
  "fewest-peptides-choice", "fast-fewest-peptides-choice"
};

static const char* hardklor_hardklor_algorithm_type_strings[NUMBER_HK_ALGORITHM_TYPES] = {
  "invalid", "Basic", "FewestPeptides", "FastFewestPeptides",
  "FewestPeptidesChoice", "FastFewestPeptidesChoice"
};

HARDKLOR_ALGORITHM_T string_to_hardklor_algorithm_type(const string& name) {
  int hk_algorithm = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING,
    hardklor_algorithm_type_strings, NUMBER_HK_ALGORITHM_TYPES);

  if (hk_algorithm < 0) {
    hk_algorithm = 0;
  }

  return (HARDKLOR_ALGORITHM_T)hk_algorithm;
}

string hardklor_hardklor_algorithm_type_to_string(
  HARDKLOR_ALGORITHM_T type
) {
  return string(hardklor_hardklor_algorithm_type_strings[type]);
}

char* ion_type_to_string(ION_TYPE_T type) {

  char* type_str = my_copy_string(ion_type_strings[type]);
  return type_str;

}


/*
 * The string version of SCORER_TYPE_T
 */
static const char* scorer_type_strings[NUMBER_SCORER_TYPES] = {
  "spscore",
  "xcorr_score",
  "evalue_score",

  "xcorr first",
  "xcorr second",

  "decoy_xcorr_qvalue",
  "decoy_xcorr_peptide_qvalue",
  "decoy_xcorr_PEP",

  "decoy_evalue_qvalue",
  "decoy_evalue_peptide_qvalue",
  "decoy_evalue_pep",
   
  "logp_weibull_xcorr",
  "logp_bonf_weibull_xcorr",
  "logp_qvalue_weibull_xcorr",
  "logp_weibull_PEP",
  "logp_peptide_qvalue_weibull",

  "percolator_score", 
  "percolator_qvalue",
  "percolator_peptide_qvalue",
  "percolator_PEP",

  "qranker_score", 
  "qranker_qvalue",
  "qranker_peptide_qvalue",
  "qranker_PEP",

  "barista_score", 
  "barista_qvalue",
  "barista_peptide_qvalue",
  "barista_PEP",

  "deltacn",
  "deltalcn",
  "by_ions_matched",
  "by_ions_total",
  "exact_pvalue",
  "refactored_xcorr",
  "res-ev score",
  "res-ev p-value",
  "combined p-value"
};

bool string_to_scorer_type(const string& name, SCORER_TYPE_T* result) {
  bool success = true;
  int scorer_type = convert_enum_type_str(name.c_str(), INVALID_ENUM_STRING, 
                                          scorer_type_strings,
                                          NUMBER_SCORER_TYPES);
  *result = (SCORER_TYPE_T)scorer_type;
  if (scorer_type < 0) {
    success = false;
  }
  return success;
}

const char* scorer_type_to_string(SCORER_TYPE_T type) {
  if ((int)type > NUMBER_SCORER_TYPES) {
    return NULL;
  }
  return scorer_type_strings[type];
}


/**
 * returns a heap allocated copy of the src string
 */
char* my_copy_string(const char* src) {
  if (src == NULL) {
    return NULL;
  }
  int length = strlen(src) +1; // +\0
  char* copy = 
    (char *)mymalloc(sizeof(char)*length);
  return strncpy(copy, src, length);  
}

/**
 * Returns copy of the src string upto the specified length.
 * Includes a null terminating character.
 * The string is heap allocated; thus, user must free.
 */
char* copy_string_part(const char* src, int length) {
  char* copy = (char*)mycalloc(length+1, sizeof(char));
  strncpy(copy, src, length);
  copy[length] = '\0';
  return copy;
}

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * compare the absolute value of the difference of two numbers with an 
 * appropriate epsilon to get relations.
 * Multiplying the epsilon by sum of the comparands adjusts the comparison 
 * to the range of the numbers, allowing a single epsilon to be used for many, 
 * or perhaps all compares.
 */
/*inline*/ int compare_float(FLOAT_T float_a, FLOAT_T float_b) {
  FLOAT_T EPSILON = PRECISION;
  FLOAT_T sum = float_a + float_b;
  if (fabsf(float_a - float_b) <= fabsf(sum)* EPSILON) { // a == b
    return 0;
  } else if ((float_a - float_b) > fabsf(sum)* EPSILON) { // a > b
    return 1;
  } else { // a < b
    return -1;
  }
}

/**
 * Compares two numbers and returns true if when rounded to the given
 * precision they are equal.  Otherwise, returns false.  
 * E.g. is_equal(0.10, 0.14, 1) -> true. is_equal(0.10, 0.15, 1) -> false
 */
bool is_equal(FLOAT_T a, FLOAT_T b, int precision) {
  a = (a * pow((FLOAT_T) 10.0, (FLOAT_T) precision)) + 0.5;
  b = (b * pow((FLOAT_T) 10.0, (FLOAT_T) precision)) + 0.5;
  return (int)a == (int)b;
}

/**
 * parses the filename and path  
 * returns an array A, with A[0] the filename and A[1] the path to the filename
 * returns A[1] NULL if only a filename was passed in
 * ex) ../../file_name => returns filename , ../../
 *     file_name => returns filename, NULL
 *\returns A heap allocated array of both filename and path
 */
char** parse_filename_path(const string& file) {
  int len = file.length();
  int end_idx = len;
  int end_path = -1;  // index of where the last "/" is located
  char* path = NULL;
  char* filename = NULL;
  char** result = (char**)mycalloc(2, sizeof(char*));

  for (; end_idx > 0; --end_idx) {
    if (file[end_idx - 1] == '/') {
      end_path = end_idx;
      break;
    }
  }
  // copy path, if there is a "/" in the file
  if (end_path != -1) {
    path = copy_string_part(file.c_str(), end_path);
  }
  // copy filename
  filename = copy_string_part(&(file.c_str()[end_idx]), len); 
  
  // set result with filename and path
  result[0] = filename;
  result[1] = path;
  
  return result;
}

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
char** parse_filename_path_extension(
  const char* file, ///< filename and path to parse -in
  const char* extension ///< extension to look for (includes leading .) --in
) {

  carp(CARP_DETAILED_DEBUG, "Given path/file %s and ext %s", file, extension);
  char** file_path_array = parse_filename_path(file);
  char* trimmed_filename = file_path_array[0];

  // look for extension
  if (extension != NULL) {
    carp(CARP_DETAILED_DEBUG, "File trimmed of path is %s", trimmed_filename);
    if (!StringUtils::EndsWith(trimmed_filename, extension)) {
      return file_path_array;  // extension not found, don't change filename
    }

    int file_len = strlen(trimmed_filename);
    int ext_len = strlen(extension);
    // after comparing the whole extension, it matched
    trimmed_filename[file_len - ext_len] = '\0';

  } else { // find the last "."
    char* dot = strrchr(trimmed_filename, '.');
    if (dot != NULL) {
      *dot = '\0';
    }
  }

  carp(CARP_DETAILED_DEBUG, "Final trimmed filename %s", trimmed_filename);


  return file_path_array;
}

/**
 * parses the filename
 * ex) ../../file_name => returns filename
 *\returns A heap allocated array of filename
 */
char* parse_filename(const char* file) {
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  // index of where the last "/" is located
  char* filename = NULL;
  
  for (; end_idx > 0; --end_idx) {
    if (strncmp(&file[end_idx - 1], "/", 1) == 0) {
      end_path = end_idx;
      break;
    }
  }
  
  // copy filename
  filename = copy_string_part(&file[end_path], len); 
  
  return filename;
}

/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* int_to_char(unsigned int i) {
  unsigned int digit = i / 10;
  char* int_string = (char*)mycalloc(digit+2, sizeof(char));
  sprintf(int_string, "%d", i);
  return int_string;
}
 
/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* signed_int_to_char(int i) {
  int digit = abs(i)/ 10;
  char* int_string = (char*)mycalloc(digit+2, sizeof(char));
  sprintf(int_string, "%d", i);
  return int_string;
}

/**
 * given two strings return a concatenated third string
 * \returns a heap allocated string that concatenates the two inputs
 */
char* cat_string(const char* string_one, const char* string_two) {
  int len_one = strlen(string_one);
  int len_two = strlen(string_two);
  
  char* result = (char*)mycalloc(len_one + len_two + 1, sizeof(char));
  strncpy(result, string_one, len_one);
  strncpy(&result[len_one], string_two, len_two);
  return result;
}

/**
 * Adds the fileroot parameter to a string as a prefix.
 */
string prefix_fileroot_to_name(const string& name) {
  string fileroot = Params::GetString("fileroot");
  return (fileroot.empty()) ? name : fileroot + '.' + name;
}

/**
 * \returns the filepath 'output_dir'/'fileroot'.'filename' 
 */
string make_file_path(
  const string& filename ///< the name of the file
  ) {
  string output_directory = Params::GetString("output-dir");
  string fileroot = Params::GetString("fileroot");

  ostringstream name_builder;
  name_builder << output_directory;

  if (output_directory[output_directory.length() - 1] != '/' ) {
    name_builder << "/";
  }

  if (!fileroot.empty()) {
    name_builder << fileroot << ".";
  }

  name_builder << filename;

  return name_builder.str();
}

/**
 * Given the path and the filename return a file with path
 * "path/filename".  Returns filename unchanged if path = NULL.
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(const char* path, const char* filename) {
  char* result = NULL;
  if (path == NULL || strlen(path) == 0) {
    result = my_copy_string(filename);
  } else {
    // TODO (BF 26-Feb-08) don't add second / if path already ends in /
    char* ready_path = cat_string(path, "/");
    result = cat_string(ready_path, filename);
    free(ready_path);
  }
  return result;
}

/**
 * returns the file size of the given filename
 */
long get_filesize(char *FileName) {
    struct stat file;
    // return file size
    if (!stat(FileName, &file)) {
      return file.st_size;
    }
    return 0;
}

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
  const string& output_folder, // Name of output folder.
  bool overwrite  // Whether or not to overwrite an existing dir 
) {
  int result = -1;
  bool path_is_directory = false;
  bool path_exists = false;
  struct stat stat_buffer;

  // Does the output directory alredy exist?
  if (stat(output_folder.c_str(), &stat_buffer)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      path_exists = false;
      path_is_directory = false;
    } else {
      // stat failed for some other reason
      carp(
        CARP_ERROR,
        "Unable to check for status of output directory '%s': %s.\n",
        output_folder.c_str(),
        strerror(errno)
      );
      result = -1;
    }
  } else {
    path_exists = true;
    path_is_directory = S_ISDIR(stat_buffer.st_mode);
  }

  if (path_exists) {
    if (!path_is_directory) {
      carp(CARP_ERROR,
        "A non-directory file named '%s' already exists,\n"
        "so that name can't be used for an output directory.\n",
        output_folder.c_str());
      result = -1;
    } else {
      if (!overwrite) {
        carp(CARP_WARNING,
          "The output directory '%s' already exists.\nExisting files will not"
          " be overwritten.", output_folder.c_str());
        result = 0;
      } else {
        carp(CARP_WARNING,
          "The output directory '%s' already exists.\nExisting files will be overwritten.",
          output_folder.c_str());
        result = 0;
      }
    }
  } else {
    // The directory doesn't exist, so we can create it.
    // Does this accomodate the case where one or more of the
    // parent directories doesn't exit?
    int dir_access = S_IRWXU + S_IRWXG + S_IRWXO;
    if (!boost::filesystem::create_directories(output_folder)) {
      // mkdir failed
      carp(CARP_ERROR, "Unable to create output directory '%s': %s.\n",
        output_folder.c_str(), strerror(errno));
      result = -1;
    } else {
      result = 0;
      carp(CARP_INFO, "Writing results to output directory '%s'.", output_folder.c_str());
    }
  }
  return result;
}

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
  vector<const char*> old_suffixes,
  const char* new_suffix,
  const char* new_path
) {
  // check the filename for the extension.  Use the first that matches
  for (size_t suffix_idx = 0; suffix_idx < old_suffixes.size(); suffix_idx++) {
    if (StringUtils::IEndsWith(filename, old_suffixes[suffix_idx])) {
      return generate_name_path(filename, old_suffixes[suffix_idx],
                                new_suffix, new_path);
    }
  }
  // if we got to here, none of the suffixes were found, so it
  // doesn't matter which we use
  return generate_name_path(filename, "", new_suffix, new_path);

}

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
) {
  carp(CARP_DEBUG, "Generate name given filename '%s', old suffix '%s', " \
       "new suffix '%s', new path '%s'", 
       filename, old_suffix, new_suffix, new_path);

  // parse path, filename, extension
  char** name_path = parse_filename_path_extension(filename, old_suffix);

  // add the extension
  char* new_name = cat_string(name_path[0], new_suffix);

  // add path to new filename
  char* full_new_name = get_full_filename(new_path, new_name);

  // cleanup
  free(name_path[0]);
  free(name_path[1]);
  free(name_path);
  free(new_name);

  carp(CARP_DEBUG, "Final name is '%s'", full_new_name);
  return full_new_name;
}
/**
 * given a fasta_file name it returns a name with the name_tag add to the end
 * Suffix may be NULL
 * format: suffix_myfasta_nameTag
 * \returns A heap allocated file name of the given fasta file
 */
char* generate_name(
  const char* fasta_filename,
  const char* name_tag,
  const char* file_extension,
  const char* suffix
) {
  int len = strlen(fasta_filename);
  int end_path = len;  // index of where the "." is located in the file
  char* name = NULL;
  char* after_suffix = NULL;
  int suffix_length = 0;
  char** file_n_path = NULL;
  int length = 0;

  // cut off the file extension if needed
  int end_idx;
  for (end_idx = len; end_idx > 0; --end_idx) {
    if (strcmp(&fasta_filename[end_idx - 1], file_extension) == 0) {
      end_path = end_idx - 1;
      break;
    }
  }
  
  // check suffix
  if (suffix != NULL) {
    suffix_length = strlen(suffix);
    file_n_path = parse_filename_path(fasta_filename);
  }

  name = (char*)mycalloc(
      suffix_length + end_path + strlen(name_tag) + 1, sizeof(char));
  after_suffix = name;
  
  // if suffix exit add to top
  if (suffix_length != 0) {
    length = strlen(file_n_path[1]);
    if (file_n_path[1] != NULL) {
      strncpy(name, file_n_path[1], length);
      after_suffix = &name[length];
    }    
    strncpy(after_suffix, suffix, suffix_length);
    after_suffix = &after_suffix[suffix_length];
    
    length = strlen(file_n_path[0]);
    
    strncpy(after_suffix, file_n_path[0], (length - (len-end_path)));
    after_suffix = &after_suffix[(length - (len-end_path))];
    
    free(file_n_path[0]);
    free(file_n_path[1]);
    free(file_n_path);
  } else {
    strncpy(after_suffix, fasta_filename, end_path);
  }
  
  strcat(after_suffix, name_tag);
  return name;
}

/**
 * checks if each AA is an AA
 *\returns true if sequence is valid else, false
 */
bool valid_peptide_sequence(const string& sequence) {
  for (string::const_iterator i = sequence.begin(); i != sequence.end(); i++) {
    if (!isupper(*i)) {
      return false;
    }
  }
  return true;
}

/**
 * \brief Open and create a file of the given name in the given
 * directory.
 *
 * Assumes the directory exists.  Fails if file can't be opened or if
 * file exists and overwrite is false.
 *\returns A file handle to the newly created file.
 */
FILE* create_file_in_path(
  const string& filename,  ///< the filename to create & open -in
  const string& directory,  ///< the directory to open the file in -in
  bool overwrite  ///< replace file (T) or die if exists (F)
) {
  char* file_full_path = get_full_filename(directory.c_str(), filename.c_str());
  // FIXME CEG consider using stat instead
  FILE* file = fopen(file_full_path, "rb"); //to test if file exists
  if (file != NULL) {
    //The file exists, are we allowed to overwrite it?
    fclose(file);
    file = NULL;
    if (!overwrite) {
        // Not allowed to overwrite, we must die.
        carp(CARP_FATAL, 
          "The file '%s' already exists and cannot be overwritten. "
          "Use --overwrite T to replace or choose a different output file name",
          file_full_path);
    } else {
      // Allowed to overwrite, send warning message.
      carp(CARP_WARNING, 
        "The file '%s' already exists and will be overwritten.", file_full_path);
    }
  }
  
  file = fopen(file_full_path, "wb+"); //read and write, replace existing

  if (file == NULL) {
    carp(CARP_FATAL, "Failed to create and open file: %s", file_full_path);
  }
  
  free(file_full_path);

  return file;
}

ofstream* create_stream_in_path(
  const char* filename,  ///< the filename to create & open -in
  const char* directory,  ///< the directory to open the file in -in
  bool overwrite  ///< replace file (T) or die if exists (F)
) {
  char* file_full_path = get_full_filename(directory, filename);
  // FIXME CEG consider using stat instead
  FILE* file = fopen(file_full_path, "rb"); //to test if file exists
  if (file != NULL) {  
    //The file exists, are we allowed to overwrite it?
    fclose(file);
    file = NULL;
    if (!overwrite) {
        // Not allowed to overwrite, we must die.
        carp(CARP_FATAL, 
          "The file '%s' already exists and cannot be overwritten. "
          "Use --overwrite T to replace or choose a different output file name",
          file_full_path);
    } else {
      // Allowed to overwrite, send warning message.
      carp(CARP_WARNING, 
        "The file '%s' already exists and will be overwritten.", file_full_path);
    }
  }
  
  ofstream* fout = new ofstream(file_full_path);

  if (fout == NULL) {
    carp(CARP_FATAL, "Failed to create and open file: %s", file_full_path);
  }
  free(file_full_path);
  return fout;
}

bool get_first_last_scan_from_string(
  const std::string& const_scans_string,
  int& first_scan,
  int& last_scan
  ) {

  set<int> scans;

  if (get_scans_from_string(const_scans_string, scans)) {
    first_scan = *(scans.begin()++);
    last_scan = *(scans.rbegin()++);
    carp(CARP_DEBUG, "scan string:%s %i %i", const_scans_string.c_str(), first_scan, last_scan);
    return true;
  } else {
    return false;
  }
}

bool get_scans_from_string(
  const string& const_scans_string,
  set<int>& scans) {
  scans.clear();

  //first tokenize by comma.
  vector<string> tokens_comma = StringUtils::Split(const_scans_string, ',');
  if (tokens_comma.size() > 1) {
    carp_once(CARP_WARNING, "Multiple scans detected in line %s. "
      "Crux currently only handles "
      "first_scan - last_scan properly", const_scans_string.c_str());
  }
  int temp_scan;

  for (size_t idx1 = 0; idx1 < tokens_comma.size(); idx1++) {
    string current = tokens_comma[idx1];
    if (current.find("-") == string::npos) {
      if (StringUtils::TryFromString(current, &temp_scan)) {
        scans.insert(temp_scan);
      } else {
        carp(CARP_ERROR, "Error parsing scans line:%s", const_scans_string.c_str());
        return false;
      }
    } else {
      vector<string> tokens_dash = StringUtils::Split(tokens_comma[idx1], '-');
      if (tokens_dash.size() != 2) {
        carp(CARP_ERROR, "Error parsing scans line:%s here:%s",
          const_scans_string.c_str(), tokens_comma[idx1].c_str());
        return false;
      }
      int temp_scan2;
      bool success = StringUtils::TryFromString(tokens_dash[0], &temp_scan);
      success &= StringUtils::TryFromString(tokens_dash[1], &temp_scan2);
      if (!success || temp_scan > temp_scan2) {
        carp(CARP_ERROR, "Error parsing scans line:%s here: %s", 
          const_scans_string.c_str(), tokens_comma[idx1].c_str());
      }
      for (int idx3 = temp_scan;idx3 <= temp_scan2 ; idx3++) {
        scans.insert(idx3);
      }
    }
  }
  return true;
}

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
) {  
  return (myrandom_limit(high - low + 1) + low);
}

/**
 *\returns the number of digits in the number
 */
int get_number_digits(
  int number ///< the number to count digits
) {
  int idx = 0;
  for (; number >= 10; ++idx) {
    number = number/10;    
  }
  return ++idx;
}

void swap_quick(
  FLOAT_T* a,
  int idx,
  int jdx
) {
  FLOAT_T temp = 0;
  temp = a[idx];
  a[idx] = a[jdx];
  a[jdx] = temp;
}
 
long Random(int i, int j) {
  return i + myrandom_limit(j-i+1);
}

void quick_sort(FLOAT_T a[], int left, int right) {
  int last = left, i;

  if (left >= right) return;
  
  swap_quick(a, left, Random(left, right));
  for (i = left + 1; i <= right; i++) {
    if (a[i] > a[left]) { /// CHECK THIS!!
      swap_quick(a, ++last, i);
    }
  }
  swap_quick(a, left, last);
  quick_sort(a, left, last-1);
  quick_sort(a, last+1, right);
}

void quicksort(FLOAT_T a[], int array_size) {
  quick_sort(a, 0, array_size - 1);
}

/**
 * \brief Shuffle an array of floats.  Uses the Knuth algorithm.  Uses
 * get_random_number_interval() to generate random numbers. 
 */
void shuffle_floats(FLOAT_T* array, int size) {
  if (array == NULL) {
    carp(CARP_ERROR, "Cannot shuffle NULL array.");
    return;
  }
  int idx, switch_idx;
  int last_element_idx = size - 1;
  FLOAT_T temp_value;
  for (idx = 0; idx < size; idx++) {
    switch_idx = get_random_number_interval(idx, last_element_idx);
    temp_value = array[idx];
    array[idx] = array[switch_idx];
    array[switch_idx] = temp_value;
  }
}

/**
 * Fits a three-parameter Weibull distribution to the input data. 
 * Implementation of Weibull distribution parameter estimation from 
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta, c (which in this case is the amount the data should
 * be shifted by) and the best correlation coefficient
 */
void fit_three_parameter_weibull(
  FLOAT_T* data, ///< the data to be fit -in
  int fit_data_points, ///< the number of data points to fit -in
  int total_data_points, ///< the total number of data points to fit -in
  FLOAT_T min_shift, ///< the minimum shift to allow -in
  FLOAT_T max_shift, ///< the maximum shift to allow -in
  FLOAT_T step,      ///< step for shift -in
  FLOAT_T corr_threshold, ///< minimum correlation, else no fit -in
  FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
  FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
  FLOAT_T* shift,     ///< the best shift -out
  FLOAT_T* correlation   ///< the best correlation -out
) {
  FLOAT_T correlation_tolerance = 0.1;
  
  FLOAT_T best_eta = 0.0;
  FLOAT_T best_beta = 0.0;
  FLOAT_T best_shift = 0.0;
  FLOAT_T best_correlation = 0.0;

  FLOAT_T cur_eta = 0.0;
  FLOAT_T cur_beta = 0.0;
  FLOAT_T cur_correlation = 0.0;
  FLOAT_T cur_shift = 0.0;

  for (cur_shift = max_shift; cur_shift > min_shift ; cur_shift -= step) {

    fit_two_parameter_weibull(data, fit_data_points, total_data_points, 
                              cur_shift, &cur_eta, &cur_beta, &cur_correlation);
    //According to the definition of the weibull distribution,
    //https://en.wikipedia.org/wiki/Weibull_distribution
    //the eta and beta parameters both have to be >0.
    //Put in a check here so that any fit needs to satisfy these constraints
    //In order to be considered.
    //SJM 2016_03_03
    if (cur_beta > 0 && cur_eta > 0 && cur_correlation > best_correlation){
      best_eta = cur_eta;
      best_beta = cur_beta;
      best_shift = cur_shift;
      best_correlation = cur_correlation;
    } else if (cur_correlation < best_correlation - correlation_tolerance) {
      break;
    }
  }

  // Only store the parameters if the fit was good enough.
  *correlation = best_correlation;
  if (best_correlation >= corr_threshold) {
    *eta = best_eta;
    *beta = best_beta;
    *shift = best_shift;
  } else {
    *eta = 0.0;
    *beta = 0.0;
    *shift = 0.0;
  }
  carp(CARP_DEBUG, "e:%g b:%g s:%g c:%g", *eta, *beta, *shift, *correlation);
}

/**
 * Fits a two-parameter Weibull distribution to the input data. 
 *
 * Called by the three parameter weibull fitting function to see if
 * the proposed shift gives the best correlation.  If there are too
 * few data points, sets correlation to 0 (minimum value).
 * http:// www.chinarel.com/onlincebook/LifeDataWeb/rank_regression_on_y.htm
 * \returns eta, beta and the correlation coefficient.
 */
void fit_two_parameter_weibull(
    FLOAT_T* data, ///< the data to be fit. should be in descending order -in
    int fit_data_points, ///< the number of data points to fit -in
    int total_data_points, ///< the total number of data points -in
    FLOAT_T shift, ///< the amount by which to shift our data -in
    FLOAT_T* eta,      ///< the eta parameter of the Weibull dist -out
    FLOAT_T* beta,      ///< the beta parameter of the Weibull dist -out
    FLOAT_T* correlation ///< the best correlation -out
) {
  FLOAT_T* X = (FLOAT_T*)mymalloc(sizeof(FLOAT_T) * fit_data_points); //hold data here

  // transform data into an array of values for fitting
  // shift (including only non-neg data values) and take log
  int idx;
  for (idx = 0; idx < fit_data_points; idx++) {
    FLOAT_T score = data[idx] + shift; // move right by shift
    if (score <= 0.0) {
      carp(CARP_DEBUG, "Reached negative score at idx %i", idx);
      fit_data_points = idx;
      break;
    } 
    X[idx] = log(score);
    // carp(CARP_DEBUG, "X[%i]=%.6f=ln(%.6f)", idx, X[idx], score);
  }

  FLOAT_T* Y   = (FLOAT_T*)mymalloc(sizeof(FLOAT_T) * fit_data_points);
  for (idx = 0; idx < fit_data_points; idx++) {
    int reverse_idx = total_data_points - idx;
    FLOAT_T F_T_idx = (reverse_idx - 0.3) / (total_data_points + 0.4);
    Y[idx] = log( -log(1.0 - F_T_idx) );
    //carp(CARP_DEBUG, "Y[%i]=%.6f", idx, Y[idx]);
  }

  int N = fit_data_points; // rename for formula's sake
  FLOAT_T sum_Y  = 0.0;
  FLOAT_T sum_X  = 0.0;
  FLOAT_T sum_XY = 0.0;
  FLOAT_T sum_XX = 0.0;
  for (idx = 0; idx < fit_data_points; idx++) {
    sum_Y  += Y[idx];
    sum_X  += X[idx];
    sum_XX += X[idx] * X[idx];
    sum_XY += X[idx] * Y[idx];
  }
  carp(CARP_DETAILED_DEBUG, "sum_X=%.6f", sum_X);
  carp(CARP_DETAILED_DEBUG, "sum_Y=%.6f", sum_Y);
  carp(CARP_DETAILED_DEBUG, "sum_XX=%.6f", sum_XX);
  carp(CARP_DETAILED_DEBUG, "sum_XY=%.6f", sum_XY);

  FLOAT_T b_num    = sum_XY - (sum_X * sum_Y / N);
  carp(CARP_DETAILED_DEBUG, "b_num=%.6f", b_num);
  FLOAT_T b_denom  = sum_XX - sum_X * sum_X / N;
  carp(CARP_DETAILED_DEBUG, "b_denom=%.6f", b_denom);
  FLOAT_T b_hat    = b_num / b_denom;

  FLOAT_T a_hat    = (sum_Y - b_hat * sum_X) / N;
  *beta = b_hat;
  *eta  = exp( - a_hat / *beta );

  FLOAT_T c_num   = 0.0;
  FLOAT_T c_denom_X = 0.0;
  FLOAT_T c_denom_Y = 0.0;
  FLOAT_T mean_X = sum_X / N;
  FLOAT_T mean_Y = sum_Y / N;
  for (idx = 0; idx < N; idx++) {
    FLOAT_T X_delta = X[idx] - mean_X; 
    FLOAT_T Y_delta = Y[idx] - mean_Y;
    c_num += X_delta * Y_delta;
    c_denom_X += X_delta * X_delta;
    c_denom_Y += Y_delta * Y_delta;
  }
  FLOAT_T c_denom = sqrt(c_denom_X * c_denom_Y);
  if (c_denom == 0.0) {
    //carp(CARP_FATAL, "Zero denominator in correlation calculation!");
    carp(CARP_DETAILED_DEBUG, "Zero denominator in correlation calculation!");
    *correlation = 0.0; // min value
    *eta = 0;
    *beta = 0;
  } else {
    *correlation = c_num / c_denom;
  }
  carp(CARP_DETAILED_DEBUG, "shift=%.6f", shift);
  carp(CARP_DETAILED_DEBUG, "eta=%.6f", *eta);
  carp(CARP_DETAILED_DEBUG, "beta=%.6f", *beta);
  carp(CARP_DETAILED_DEBUG, "correlation=%.6f", *correlation);

  free(Y);
  free(X);
}


/**
 * \brief Open the fasta file and prepare it for
 * searching.  Die if the input file cannot be found or read.
 * \returns The number of proteins in the file
 */
int prepare_protein_input(
  const string& input_file,     ///< name of the fasta file
  Database** database ///< return new fasta database here
) {
  int num_proteins = 0;

  carp(CARP_INFO, "Preparing protein fasta file %s", input_file.c_str());
  *database = new Database(input_file, false);         
  if (database == NULL) {
    carp(CARP_FATAL, "Could not create protein database");
  } 

  if (!(*database)->parse()) {
    carp(CARP_FATAL, "Error with protein database.");
  } 
  return (*database)->getNumProteins();
}

/**
 * Maximum characters per line when printing formatted text.
 */
static const int MAX_CHARS_PER_LINE = 70;

/**
 *\brief Extend a given string with lines not exceeding a specified width, 
 * breaking on spaces.
 */
void strcat_formatted(
  char*       string_to_extend,
  const char* lead_string,        // Appears at the start of each line.
  const char* extension           // Text to add.
) {
  int i_string = 0;
  int string_length = strlen(extension);
  char buffer[MAX_CHARS_PER_LINE + 1];

  while (string_length - i_string > MAX_CHARS_PER_LINE) {

    // Be sure to break on a space.
    int index_of_last_space = MAX_CHARS_PER_LINE;
    while (extension[i_string + index_of_last_space] != ' ') {
      index_of_last_space--;
    }

    sprintf(buffer, "%.*s", index_of_last_space, &(extension[i_string]));

    strcat(string_to_extend, lead_string);
    strcat(string_to_extend, buffer);
    strcat(string_to_extend, "\n");

    i_string += index_of_last_space + 1;
  }
  strcat(string_to_extend, lead_string);
  strcat(string_to_extend, &(extension[i_string]));
  strcat(string_to_extend, "\n");
}

/**
 * \brief Checks if the given input file contains target, decoy PSMs or 
 * concatenated search results.
 *
 *\returns corrected file names. It does not check if files are exist.
 */
void check_target_decoy_files(
  string &target,   //filename of the target PSMs
  string &decoy     //filename of the decoy PSMs
) {
  int target_pos = target.find("target");
  if (target_pos < 0) {
    int decoy_pos = decoy.find("decoy");
    if (decoy_pos < 0) {
      // user gave concatenated result file
      decoy = "";
    } else {
      // user gave decoy results file
      target.replace(decoy_pos, DECOY_STRING_LENGTH, "target");
    }
  } else {
    // user gave target results file
    decoy.replace(target_pos, TARGET_STRING_LENGTH, "decoy"); 
  }
}

void get_search_result_paths(
  const vector<string>& infiles,
  vector<string>& outpaths ///< paths of all search results -out 
) {
  outpaths.clear();
  for (vector<string>::const_iterator i = infiles.begin(); i != infiles.end(); i++) {
    if (!FileUtils::Exists(*i)) {
      carp(CARP_ERROR, "Search file '%s' doesn't exist", i->c_str());
      continue;
    }
    outpaths.push_back(*i);
  }
  for (vector<string>::iterator i = outpaths.begin(); i != outpaths.end(); i++) {
    string path = *i;
    string target = path;
    string decoy = path;
    check_target_decoy_files(target, decoy);
    if (!target.empty() && target != path && FileUtils::Exists(target) &&
        find(outpaths.begin(), outpaths.end(), target) == outpaths.end()) {
      i = outpaths.insert(i + 1, target);
    }
    if (!decoy.empty() && decoy != path && FileUtils::Exists(decoy) &&
        find(outpaths.begin(), outpaths.end(), decoy) == outpaths.end()) {
      i = outpaths.insert(i + 1, decoy);
    }
  }
}

void get_files_from_list(
  const string &infile, ///< path of the first file.
  std::vector<std::string> &outpaths ///< paths of all search results -out
  ) {
  outpaths.clear();
  if (Params::GetBool("list-of-files")) {
    LineFileReader reader(infile);
    while (reader.hasNext()) {
      string current = reader.next();
      carp(CARP_INFO, "current is:%s", current.c_str());
      if (FileUtils::Exists(current)) {
        outpaths.push_back(current);
      } else {
        carp(CARP_ERROR, "Search file '%s' doesn't exist", current.c_str());
      }
    }
  } else if (FileUtils::Exists(infile)) {
    outpaths.push_back(infile);
  }
}

bool parseUrl(string url, string* host, string* path) {
  if (!host || !path) {
    return false;
  }
  // find protocol
  size_t protocolSuffix = url.find("://");
  if (protocolSuffix != string::npos) {
    url = url.substr(protocolSuffix + 3);
  }
  size_t pathBegin = url.find('/');
  if (pathBegin == string::npos) {
    *host = url;
    *path = "/";
  } else {
    *host = url.substr(0, pathBegin);
    *path = url.substr(pathBegin);
  }
  if (host->empty()) {
    *host = *path = "";
    return false;
  }
  return true;
}

string httpRequest(const string& url, const std::string& data, bool waitForResponse) {
  // Parse URL into host and path components
  string host, path;
  if (!parseUrl(url, &host, &path)) {
    carp(CARP_ERROR, "Failed parsing URL %s", url.c_str());
    return "";
  }

  using namespace boost::asio;

  // Establish TCP connection to host on port 80
  io_service service;
  ip::tcp::resolver resolver(service);
  ip::tcp::resolver::iterator endpoint = resolver.resolve(ip::tcp::resolver::query(host, "80"));
  ip::tcp::socket sock(service);
  connect(sock, endpoint);

  // Determine method (GET if no data; otherwise POST)
  string method = data.empty() ? "GET" : "POST";
  string contentLengthHeader = data.empty()
    ? ""
    : "Content-Length: " + StringUtils::ToString(data.length()) + "\r\n";
  // Send the HTTP request
  string request =
    method + " " + path + " HTTP/1.1\r\n"
    "Host: " + host + "\r\n" +
    contentLengthHeader +
    "Connection: close\r\n"
    "\r\n" + data;
  sock.send(buffer(request));

  if (!waitForResponse) {
    return "";
  }

  // Read the response into stringstream
  stringstream response;
  boost::system::error_code ec;
  do {
    char buf[1024];
    size_t n = sock.receive(buffer(buf), 0, ec);
    if (!ec) {
      response.write(buf, n);
    }
  } while (!ec);
  return response.str();
}

void postToAnalytics(const string& appName) {
  // Post data to Google Analytics
  // For more information, see: https://developers.google.com/analytics/devguides/collection/protocol/v1/devguide
  try {
    stringstream paramBuilder;
    paramBuilder << "v=1"                // Protocol verison
                 << "&tid=UA-26136956-1" // Tracking ID
                 // TODO Generate UUID for cid
                 // The Client ID (cid) anonymously identifies a particular user
                 // or device and should be a random UUID
                 << "&cid=35009a79-1a05-49d7-b876-2b884d0f825b"
                 << "&t=event"           // Hit type
                 << "&ec=crux"           // Event category
                 << "&ea=" << appName    // Event action
                 << "&el="               // Event label
#ifdef _MSC_VER
                      "win"
#elif __APPLE__
                      "mac"
#else
                      "linux"
#endif
                   << '-' << CRUX_VERSION;
      httpRequest(
        "http://www.google-analytics.com/collect",
        paramBuilder.str(),
        false);
  } catch (...) {
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

