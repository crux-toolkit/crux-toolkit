/******************************************************************************
 * FILE: parameter-file.h
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities.
 *****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include "utils.h"

#define PARAMETER_LENGTH 1024
#define NUM_PARAMS 512
#define MAX_LINE_LENGTH 4096

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  );

/**
 * 
 * Each of the following functions searches through the list of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value, or the default value if
 * the parameter is not found.
 * \returns TRUE if paramater value or default is TRUE, else FALSE
 */ 
BOOLEAN_T get_boolean_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T default_value  ///< the dafault value to use if not found -in
 );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter_array, and the dafault otherwise.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter or default value
 */
int get_int_parameter(
  char* name,  ///< the name of the parameter looking for -in
  int default_value  ///< the dafault value to use if not found -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter_array, and the dafault otherwise.  This
 * function exits if there is a conversion error. 
 *\returns the double value of the parameter of the default value
 */
double get_double_parameter(
  char* name,   ///< the name of the parameter looking for -in
  double default_value   ///< the dafault value to use if not found -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, return NULL.
 * \returns the string value to which matches the parameter name, else returns NULL
 */
char* get_string_parameter(
  char* name  ///< the name of the parameter looking for -in
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

/**
 * Check to see if any parameters were not used.  Issue a warning.
 */
void check_unused_parameters(void);

#endif
