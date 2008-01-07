/**
 * \file parameter.h
 * $Revision: 1.26 $
 * \brief General parameter handling utilities. MUST declare ALL optional command line parameters here inside initalialize_parameters
 ****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "utils.h"
#include "crux-utils.h"
#include "carp.h"
#include "hash.h"
#include "objects.h"
#include "peptide.h"
#include "spectrum.h"
#include "peak.h"
#include "mass.h"
#include "scorer.h"
#include "parse_arguments.h"

#define PARAMETER_LENGTH 1024 ///< default length of parameter name and value in characters
#define NUM_PARAMS 512 ///< initial number of parameters allowed
#define MAX_LINE_LENGTH 4096 ///< maximum length of a line on the parameter file
#define BILLION 1000000000.0
#define SMALL_BUFFER 256
#define MAX_SET_PARAMS 256

#define NUMBER_PARAMETER_TYPES 11
enum parameter_type {INT_P, DOUBLE_P, STRING_P, MASS_TYPE_P, 
		     PEPTIDE_TYPE_P, BOOLEAN_P, SORT_TYPE_P,
		     SCORER_TYPE_P, OUTPUT_TYPE_P, ION_TYPE_P, ALGORITHM_TYPE_P};
typedef enum parameter_type PARAMETER_TYPE_T;

/**
 * initialize parameters
 * Every required argument for every executable 
 * and every option and its default value
 * must be declared here
 */
void initialize_parameters(void);

/**
 * free heap allocated parameters hash table
 */
void free_parameters(void);

/**
 * Identify which of the parameters can be changed 
 * on the command line.  Provide a list of the parameter names
 * and the number of parameters in that list.
 * Requires that initialize_parameters() has been run.
 */
BOOLEAN_T select_cmd_line_options(char**, int);

/**
 * Identify what are the required arguments
 * on the command line.  Provide a list of the argument names
 * and the number of arguments in that list.
 * Requires that initialize_parameters() has been run.
 */
BOOLEAN_T select_cmd_line_arguments(char**, int);

/**
 * Take the command line string from main, find the parameter fil
 * option (if present), parse it's values into the hash, and parse
 * the command line options and arguments into the hash
 * main then retrieves the values through get_value
 */
BOOLEAN_T parse_cmd_line_into_params_hash(int, char**, char*);

/**
 * Each of the following functions searches through the hash table of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value.
 * \returns TRUE if paramater value is TRUE, else FALSE
 */ 
BOOLEAN_T get_boolean_parameter(
 char*     name  ///< the name of the parameter looking for -in
 );

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter
 */
int get_int_parameter(
  char* name  ///< the name of the parameter looking for -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the double value of the parameter
 */
double get_double_parameter(
  char* name   ///< the name of the parameter looking for -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, abort.
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter(
  char* name  ///< the name of the parameter looking for -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is a pointer to the original string
 * Thus, user should not free, good for printing
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter_pointer(
  char* name  ///< the name of the parameter looking for -in
  );

MASS_TYPE_T get_mass_type_parameter(
 char* name
 );

SORT_TYPE_T get_sort_type_parameter(
 char* name
 );

ALGORITHM_TYPE_T get_algorithm_type_parameter(
 char* name
 );

SCORER_TYPE_T get_scorer_type_parameter(
 char* name
 );

MATCH_SEARCH_OUTPUT_MODE_T get_output_type_parameter(
 char* name
 );

ION_TYPE_T get_ion_type_parameter(
 char* name
 );

/**
 * Searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * Returns a peptide_type enumerated type (in objects.h)
 */ 
PEPTIDE_TYPE_T get_peptide_type_parameter(
  char* name
  );

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


#endif
