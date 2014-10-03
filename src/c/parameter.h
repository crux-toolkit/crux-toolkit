/**
 * \file parameter.h
 * $Revision: 1.37 $
 * \brief General parameter handling utilities. All values stored here.

 * \detail MUST declare ALL optional command line parameters and
 * required command line arguments here in initalialize_parameters.
 * Parameters can only be set via parse_cmd_line_into_params_hash()
 * which takes values from the command line and from an optional
 * parameter file (provided on the command line with the --parameter
 * option). Options are checked for correct type and legal values.
 * Exits with usage statement on error.  Parameter values can be
 * retrieved with get_<type>_paramter functions.
 * 
 ****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <string>
#include "utils.h"
#include "crux-utils.h"
#include "carp.h"
#include "hash.h"
#include "objects.h"

#include "parse_arguments.h"
#include "modifications.h"

static const int PARAMETER_LENGTH = 1024; 
///< maximum length of parameter name and value in characters
static const int NUM_PARAMS = 512; ///< maximum number of parameters allowed
static const int MAX_LINE_LENGTH = 4096; ///< maximum line length in the parameter file
static const int BILLION = 1000000000;
static const int SMALL_BUFFER = 256;
static const int MAX_SET_PARAMS = 256;

// The size of the bins for discretizing the m/z axis of the
// observed spectrum.  For use with monoisotopic mass.
static const FLOAT_T BIN_WIDTH_MONO = 1.0005079;
// The size of the bins for discretizing the m/z axis of the
// observed spectrum.  For use with average mass.
static const FLOAT_T BIN_WIDTH_AVERAGE = 1.0011413;

// Global variables
// NOTE (BF mar-10-09): Could be like mod lists, but will require a
// get_parameter call for each...residue? protein?
extern char* pre_cleavage_list;
extern char* post_cleavage_list;
extern int pre_list_size;
extern int post_list_size;
extern bool pre_for_inclusion;
extern bool post_for_inclusion;

// TODO (BF 1-28-08): these should be private. move to parameter.c
/**
 * Data types of parameters.  Used for checking valid parameter input
 * from user.
 *
 * To add a new parameter type:  
 *  (for NEW types) 1. create enum, 
 *                  2. create array of strings, 
 *                  3. write string-to-type,
 *                  4. write type-to-string
 *  (for ALL types) 5. add to parameter-type enum and strings
 *                  6. write get-type-parameter
 *                  7. write set-type-parameter
 *                  8. add to the check_type_and_bounds
 *
 */
enum parameter_type {
  INT_P,             ///< parameters of type int
  DOUBLE_P,          ///< parameters of type double
  STRING_P,          ///< parameters of type char*
  MASS_TYPE_P,       ///< parameters of type MASS_TYPE_T
  DIGEST_TYPE_P,   ///< parameters of type DIGEST_T
  ENZYME_TYPE_P,     ///< parameters of type ENZYME_T
  //PEPTIDE_TYPE_P,    ///< parameters of type PEPTIDE_TYPE_T
  BOOLEAN_P,         ///< parameters of type bool
  SCORER_TYPE_P,     ///< parameters of type SCORER_TYPE_T
  ION_TYPE_P,        ///< parameters of type ION_TYPE_T
  ALGORITHM_TYPE_P,  ///< parameters of type ALGORITHM_TYPE_T
  HARDKLOR_ALGORITHM_TYPE_P, ///< parameters of type HARDKLOR_ALGORITHM_T
  SPECTRUM_PARSER_P, ///< parameters of type SPECTRUM_PARSER_T
  WINDOW_TYPE_P,     ///< parameters of type WINDOW_TYPE_T
  MEASURE_TYPE_P,    ///< parameters of type MEASURE_TYPE_T
  THRESHOLD_P,       ///< parameters of type THRESHOLD_TYPE_T
  PARSIMONY_TYPE_P,  ///< parameters of type PARSIMONY_TYPE_T
  QUANT_LEVEL_TYPE_P,///< parameters of type QUANT_LEVEL_TYPE_T
  DECOY_TYPE_P,      ///< parameters of type DECOY_TYPE_T
  MASS_FORMAT_P,     ///< parameters of type MASS_FORMAT_T

  NUMBER_PARAMETER_TYPES  ///< leave this last, number of types
};
typedef enum parameter_type PARAMETER_TYPE_T;

/**
 * /brief Initialize parameters to default values.
 *
 * Every required argument for every executable and every option and
 * its default value must be declared here.  Allocates hash tables for
 * holding parameter values.
 */
void initialize_parameters(void);

/**
 * free heap allocated parameters hash table
 */
void free_parameters(void);

/**
 * /brief Identify which of the parameters can be changed 
 * on the command line.  
 * /detail Provide a list of the parameter names
 * and the number of parameters in that list.
 * Requires that initialize_parameters() has been run.
 * /returns TRUE on success.
 */
bool select_cmd_line_options(const char**, int);

/**
 * /brief Identify the required command line arguments.
 * /detail  Provide a list of the argument names
 * and the number of arguments in that list.
 * Requires that initialize_parameters() has been run.
 * /returns TRUE on success.
 */
bool select_cmd_line_arguments(const char**, int);

/**
 * Take the command line string from main, find the parameter fil
 * option (if present), parse it's values into the hash, and parse
 * the command line options and arguments into the hash
 * main then retrieves the values through get_value
 */
bool parse_cmd_line_into_params_hash(int, char**, const char*);

/**
 * Each of the following functions searches through the hash table of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value.
 * \returns TRUE if paramater value is TRUE, else FALSE
 */ 
bool get_boolean_parameter(
 const char*     name  ///< the name of the parameter looking for -in
 );

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter
 */
int get_int_parameter(
  const char* name  ///< the name of the parameter looking for -in
  );

std::vector<int> get_int_vector_parameter(
  const char* name ///< the name of the parameter looking for -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the double value of the parameter
 */
double get_double_parameter(
  const char* name   ///< the name of the parameter looking for -in
  );

std::vector<double> get_double_vector_parameter(
  const char* name
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, abort.
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter(
  const char* name  ///< the name of the parameter looking for -in
  );

std::vector<std::string> get_string_vector_parameter(
  const char* name
  );
/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is a pointer to the original string
 * Thus, user should not free, good for printing
 * \returns the string value to which matches the parameter name, else aborts
 */
const char* get_string_parameter_pointer(
  const char* name  ///< the name of the parameter looking for -in
  );

MASS_TYPE_T get_mass_type_parameter(
 const char* name
 );

char get_delimiter_parameter(
  const char* name
  );

ALGORITHM_TYPE_T get_algorithm_type_parameter(
 const char* name
 );

SCORER_TYPE_T get_scorer_type_parameter(
 const char* name
 );

ION_TYPE_T get_ion_type_parameter(
 const char* name
 );

DIGEST_T get_digest_type_parameter(
  const char* name
  );

ENZYME_T get_enzyme_type_parameter(
  const char* name
  );

DIGEST_T get_digest_type_parameter(
  const char* name
  );

ENZYME_T get_enzyme_type_parameter(
  const char* name
  );

WINDOW_TYPE_T get_window_type_parameter(
  const char* name
  );

THRESHOLD_T get_threshold_type_parameter(
  const char* name
  );

PARSIMONY_TYPE_T get_parsimony_type_parameter(
  const char* name
  );

QUANT_LEVEL_TYPE_T get_quant_level_type_parameter(
  const char* name
  );

MEASURE_TYPE_T get_measure_type_parameter(
  const char* name
  );

DECOY_TYPE_T get_decoy_type_parameter(
  const char* name
  );

DECOY_TYPE_T get_tide_decoy_type_parameter(
  const char* name
  );

MASS_FORMAT_T get_mass_format_type_parameter(
  const char* name
  );

int get_max_ion_charge_parameter(
  const char* name
  );

HARDKLOR_ALGORITHM_T get_hardklor_algorithm(
  const char* name
  );

SPECTRUM_PARSER_T get_spectrum_parser_parameter(
  const char* name
  );


double get_mz_bin_width();
 
double get_mz_bin_offset();

COLTYPE_T get_column_type_parameter(
  const char* name
  );

COMPARISON_T get_comparison_parameter(
  const char* name
  );

/**
 * \returns the comet enzyme info lines parsed from the file
 * or generated defaults
 */
const std::vector<std::string>& get_comet_enzyme_info_lines();

/**
 * \brief prints all parameters except mods into the output stream
 * in xml format. 
 *
 * Each parameter has a self closing tag and has attributes name for 
 * parameter name and value for parameter value
 */
void print_parameters_xml(FILE* output);


/**
 * Prints the parameters.  If lead_string is not null, preprends it to
 * each line.
 */
void print_parameters(
  char* first_line,  ///< the first line to be printed before the parameter list -in
  char* parameter_filename,  ///< the parameter file name -in
  char* lead_string,  ///< the lead string to be printed before each line -in
  FILE* outstream  ///< the output stream -out
  );


/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Does not include the c- and n-term mods.  Argument is a
 * reference to an array of pointers.  Return 0 and set mods == NULL
 * if there are no aa_mods.
 * \returns The number of items pointed to by mods
 */
int get_aa_mod_list(AA_MOD_T*** mods);

/**
 * Adds 1 to the number of variable modifications in the
 * parameters
 */
void incrementNumMods();

/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * c-terminus.  Argument is a reference to an array of
 * pointers. Return 0 and set mods==NULL if there are no c-term 
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_c_mod_list(AA_MOD_T*** mods);

/**
 * \brief Get the pointer to the list of AA_MODs for the peptide
 * n-terminus.  Argument is a reference to an array of
 * pointers. Return 0 and set mods==NULL if there are no n-term 
 * mods.
 *
 * \returns The number of items pointed to by mods
 */
int get_n_mod_list(AA_MOD_T*** mods);

/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.  Includes aa_mods, c- and n-term mods.  Argument is a
 * reference to an array of pointers.  Returns 0 and sets mods == NULL
 * if there are no aa_mods.
 * \returns The number of items pointed to by mods
 */
int get_all_aa_mod_list(AA_MOD_T*** mods);

/**
 * \returns The index of the C_TERM or N_TERM fixed modification in
 * the global list of modifications.
 */
int get_fixed_mod_index(MOD_POSITION_T p);

/**
 * \returns the number of fixed terminal modifications: 0, 1, or 2.
 */
int get_num_fixed_mods();

/**
 * \brief Creates a file containing all parameters and their current
 * values in the parameter file format. Created in the output directory
 * named by the parameter "output-dir".
 */
void print_parameter_file(char** filename);

#endif
