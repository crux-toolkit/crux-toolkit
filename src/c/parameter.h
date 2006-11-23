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
 * Parameter file must be parsed first!
 * searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * The function sets the corresponding value,
 * if the parameter is not found or parameters has already been confirmed don't change
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_boolean_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T set_value  ///< the value to be set -in
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
 * Parameter file must be parsed first!
 * searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * The function sets the corresponding value,
 * if the parameter is not found or parameters has already been confirmed don't change
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_int_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 int set_value  ///< the value to be set -in
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
 * Parameter file must be parsed first!
 * searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * The function sets the corresponding value,
 * if the parameter is not found or parameters has already been confirmed don't change
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_double_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 double set_value  ///< the value to be set -in
 );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, abort.
 * \returns the string value to which matches the parameter name, else abort.
 */
char* get_string_parameter(
  char* name  ///< the name of the parameter looking for -in
  );

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is a pointer to the original string
 * Thus, user should no free, good for printing
 * If the value is not found, abort.
 * \returns the string value to which matches the parameter name, else abort.
 */
char* get_string_parameter_pointer(
  char* name  ///< the name of the parameter looking for -in
  );

/**
 * Parameter file must be parsed first!
 * searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * The function sets the corresponding value,
 * if the parameter is not found or parameters has already been confirmed don't change
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_string_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 char* set_value  ///< the value to be set -in
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


/**
 * set the parameters to confirmed, thus blocks any additional changes to the parameters
 */
void parameters_confirmed(void);

/**
 * Parameter file must be parsed first!
 * searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * then the function sets the corresponding value.
 * if the parameter is not found, return FALSE
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T set_options_command_line(
  char*     name,  ///< the name of the parameter looking for -in
  char* set_value, ///< the value to be set -in
  BOOLEAN_T required ///< is this a required option -in
  );

/**
 * This method should be called only after parsed command line
 * first, parse paramter file
 * Next, updates the parameter files with command line options
 * command line arguments have higher precedence 
 * parse the parameter file given the filename
 */
void parse_update_parameters(
  char* parameter_file ///< the parameter file to be parsed -in
  );

#endif
