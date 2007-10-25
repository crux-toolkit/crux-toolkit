/**
 * \file parameter.h
 * $Revision: 1.14 $
 * \brief General parameter handling utilities. MUST declare ALL optional command line parameters here inside initalialize_parameters
 *****************************************************************************/
#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H
#include "utils.h"
#include "objects.h"

#define PARAMETER_LENGTH 1024 ///< default length of parameter name and value in characters
#define NUM_PARAMS 512 ///< initial number of parameters allowed
#define MAX_LINE_LENGTH 4096 ///< maximum length of a line on the parameter file
#define BILLION 1000000000.0

/**
 * free heap allocated parameters hash table
 */
void free_parameters(void);

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
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter
 */
int get_int_parameter(
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
BOOLEAN_T set_int_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 int set_value  ///< the value to be set -in
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
 * searches through the hash table of parameters, 
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

/**
 * add parameters to parameter list
 * Should use this only for testing cases.
 * Otherwiase use normal parameter system.
 */
BOOLEAN_T add_parameter(
  char*     name,  ///< the name of the parameter to add -in
  char* set_value  ///< the value to be added -in                  
  );

/**
 * Parameter file must be parsed first!
 * Searches through the list of parameters, 
 * looking for one whose name matches the string.  
 * Returns a peptide_type enumerated type (in objects.h)
 */ 
PEPTIDE_TYPE_T get_peptide_type_parameter(
  char* name
  );

#endif
