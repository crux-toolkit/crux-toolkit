/******************************************************************************
 * FILE: parameter.c
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities. MUST declare ALL optional command parameters here inside initalialize_parameters
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "objects.h"
#include "spectrum.h"
#include "peak.h"
#include "carp.h"
#include "mass.h"
#include "scorer.h"
#include "utils.h"
#include "parameter.h"
#include "parse_arguments.h"

/**
 *\struct parameter
 *\brief the structure that handles the one parameter value
 * The parameters are stored as a global variable.
 */
struct parameter{
  BOOLEAN_T used;  ///< has this parameter been used, in other words has the user extracted the information
  char parameter_name[PARAMETER_LENGTH];  ///< the name of the parameter
  char parameter_value[PARAMETER_LENGTH]; ///< the actual parameter value corresponding to the name
};

/**
 *\struct parameter_array
 *\brief the array that holds all the different parameters
 */
struct parameter_array{
  int num_parameters;   ///< number of the total number of parameters
  struct parameter parameters[NUM_PARAMS]; ///< the paraters
};

/**
 * Global variable
 */

struct parameter_array parameters;
BOOLEAN_T parameter_initialized = FALSE; //have the parameters been initialized?
BOOLEAN_T parameter_parsed = FALSE; //have I parsed the parameter file?
BOOLEAN_T parameter_plasticity = TRUE; //can the parameters be changed?

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  );

/**
 * initialize parameters
 * ONLY add optional parameters here!!!
 * MUST declare ALL optional parameters in array to be used
 */
void initialize_parameters(void){

  //check if parameters been initialized
  if(parameter_initialized){
    carp(CARP_ERROR, "parameters has already been initialized");
    return;
  }
  
  //set number of parameters to zero
  parameters.num_parameters = 0;

  //set verbosity
  set_int_parameter("verbosity", CARP_ERROR);

  //set perameters
  set_string_parameter("parameter-file", "crux_parameter");

  //generate_peptide perameters
  set_double_parameter("min-mass", 200);
  set_double_parameter("max-mass", 2400);
  set_int_parameter("min-length", 6);
  set_int_parameter("max-length", 50);
  set_string_parameter("cleavages", "tryptic");
  set_string_parameter("isotopic-mass","average");
  set_string_parameter("redundancy", "redundant");
  set_string_parameter("use-index", "F");
  set_string_parameter("sort", "none");      // mass, length, lexical, none  
  //set_string_parameter("fasta-file", "NULL");
  set_boolean_parameter("output-sequence", FALSE);
  set_boolean_parameter("missed-cleavages", FALSE);
  
  //score_peptide_spectrum perameters
  set_double_parameter("beta", 0.075);
  set_double_parameter("max-mz", 4000);
  set_int_parameter("charge", 2);
  set_string_parameter("score-type", "sp"); 

  //now we have initialized the parameters
  parameter_initialized = TRUE;
}

/********************************************************************
 *
 * the parameter file is assumed to consist of name/value pairs,
 * separated by an equals sign.
 *
 *******************************************************************/

/**
 * add parameters to parameter list
 */
BOOLEAN_T add_parameter(
  char*     name,  ///< the name of the parameter to add -in
  char* set_value  ///< the value to be added -in                  
  )
{
  //but is there space?
  if(parameters.num_parameters >= NUM_PARAMS){
    carp(CARP_ERROR, "no more space for any additional paramerters");
    return FALSE;
  }

  // copy the name/value pairs to the right parameter
  parameters.parameters[parameters.num_parameters].used = FALSE;
  strncpy(parameters.parameters[parameters.num_parameters].parameter_name, 
          name,
          PARAMETER_LENGTH);
  strncpy(parameters.parameters[parameters.num_parameters].parameter_value,
          set_value,
          PARAMETER_LENGTH);
  parameters.num_parameters++;

  return TRUE;
}

/**
 * copy parameters to parameter list
 * there must be a matching name in the parameter list
 */
BOOLEAN_T copy_parameter(
  char*     name,  ///< the name of the parameter to add -in
  char* set_value  ///< the value to be added -in                  
  )
{
  
  int idx;
  
  //check if parameters has been initlialized
  if(!parameter_initialized){
    carp(CARP_ERROR, "must inilialize parameters before copying");
    return FALSE;
  }
    
  //check if parameters can be changed
  if(!parameter_plasticity){
   carp(CARP_ERROR, "can't change parameters once they are confirmed");
   return FALSE;
  }
  
  //check if parameter name already exist?
  for(idx = 0; idx < parameters.num_parameters; idx++){
    //if exist ovewrite it!
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      strcpy(parameters.parameters[idx].parameter_value, set_value);
      return TRUE;
    }	
  }
  
  carp(CARP_ERROR, "incorrect parameter name: %s, from the parameter file", name);
  return FALSE;
}

/**
 * This method should be called only after parsed command line
 * first, parse paramter file
 * Next, updates the parameter files with command line options
 * command line arguments have higher precedence 
 * parse the parameter file given the filename
 */
void parse_update_parameters(
  char* parameter_file ///< the parameter file to be parsed -in
  )
{
  //initialize
  initialize_parameters();
  
  //if no parameter file name has been specified in command line,
  // use default parameter filename
  if(parameter_file != NULL){
    //parse parameter file
    parse_parameter_file(parameter_file);
  }
  else{
    //parse parameter file
    parse_parameter_file(get_string_parameter_pointer("parameter-file"));
  }
  
  //update the parameters if any comman line arguments exist
  if(!update_parameter()){
    fprintf(stderr, "failed to combine command line arguemnts and parameter file\n");
    exit(-1);
  }
}

/**
 *
 * parse the parameter file given the filename
 */
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  )
{
  FILE *f;
  char *line;
  int idx;

  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    exit(-1);
  }

  //check if parameter file exist, if not exit use default parameters
  if(access(parameter_filename, F_OK)){
     carp(CARP_INFO, "no parameter_file, using default parameters");
     return;
  }

  line = (char*)mycalloc(MAX_LINE_LENGTH, sizeof(char));


  f = fopen(parameter_filename, "r");
  if(f == NULL){
    printf("couldn't open file: %s\n", parameter_filename);
    exit(-1);
  }

  while(fgets(line, MAX_LINE_LENGTH, f)==line){
    // if a line begins with "#", it is ignored as a comment 
    if(line[0] != '#'){
    /* find the '=' in the line.  Exit with error if the line 
       has no equals sign. */
      idx = 0;

      // Change the newline to a '\0'
      for(idx = MAX_LINE_LENGTH - 1; idx > 0; idx--){
	if(line[idx] == '\n' || line[idx] == '\r' || line[idx] == '\f')
	  line[idx] = '\0';
      }
      while(idx < (int)strlen(line) && line[idx] != '='){
	idx++;
      }
      if(idx == 0 || idx >= (int)(strlen(line)-1)){
	printf("lines in a parameter file must have the form:\n");
	printf("name=value\n");
	printf("in file %s, the line:\n%s\ndoes not have this format\n",
	       parameter_filename, line);
	exit(-1);
      }
      line[idx] = '\0';
      
      // copy the name/value pairs to the right parameter
      if(!copy_parameter(line, &(line[idx+1]))){
        exit(-1);
      }
    }
  }

  // Tell the user what we found.
  /*
  if(get_verbosity_level() >= CARP_INFO) {
    print_parameters("Using parameters:", parameter_filename, "  ", stderr);
  }
  */

  fclose(f);
  myfree(line);

  //now we have parsed the parameter file
  parameter_parsed = TRUE;
}

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
 )
{
  int idx;
  static char buffer[PARAMETER_LENGTH];

  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet, using default value");
    return(default_value);
  }

  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      // make sure that there is enough storage allocated in the string
      if((int)strlen(parameters.parameters[idx].parameter_value) 
	 > PARAMETER_LENGTH) {
	die("parameter %s with value %s was too long to copy to string\n",
	    parameters.parameters[idx].parameter_name,
	    parameters.parameters[idx].parameter_value);
      }
      strncpy(buffer,
	      parameters.parameters[idx].parameter_value,
	      PARAMETER_LENGTH);
      if (strcmp(buffer, "TRUE") == 0) {
	return(TRUE);
      } else if (strcmp(buffer, "FALSE") == 0) {
	return(FALSE);
      } else {
	die("Invalid Boolean parameter %s.\n", buffer);
      }
    }
  }
  carp(CARP_ERROR, "parameter name: %s, doesn't exit, using default value: %d (1=TRUE, 0=FALSE)", name, default_value);
  return(default_value);
}

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
 )
{
  int idx;
  BOOLEAN_T result;
    
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }

  //only check if parameter file has already been parsed
  if(parameter_parsed){
    //check if parameter name already exist?
    for(idx = 0; idx < parameters.num_parameters; idx++){
      //if exist ovewrite it!
      if(!strcmp(parameters.parameters[idx].parameter_name, name)){
        //set to TRUE
        if(set_value){
          strcpy(parameters.parameters[idx].parameter_value, "TRUE");
        }
        //set to FALSE
        else{
          strcpy(parameters.parameters[idx].parameter_value, "FALSE");
        }
        return TRUE;
      }	
    }
  }

  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  if(set_value){
    result = add_parameter(name, "TRUE");
  }
  else{
    result = add_parameter(name, "FALSE");
  }
  
  return result;
}

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
  )
{
  int idx;
  char *endptr;
  long int value;

  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet, using default value");
    return(default_value);
  }

  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){

      /* there is a parameter with the right name.  Now 
	 try to convert it to a base 10 integer*/
      value = strtol(parameters.parameters[idx].parameter_value, &endptr, 10);
      if ((value == LONG_MIN) || 
	  (value == LONG_MAX) || 
	  (endptr == parameters.parameters[idx].parameter_value)) {
	die("Conversion error when trying to convert parameter %s with value %s to an int\n",
	       parameters.parameters[idx].parameter_name, 
	       parameters.parameters[idx].parameter_value);
      } else {
	parameters.parameters[idx].used = TRUE;
	return((int)value);
      }
    }
  }
  carp(CARP_ERROR, "parameter name: %s, doesn't exit, using default value: %d", name, default_value);
  return(default_value);
}

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
 )
{
  int idx;
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //only check if parameter file has already been parsed
  if(parameter_parsed){
    //check if parameter name already exist?
    for(idx = 0; idx < parameters.num_parameters; idx++){
      //if exist ovewrite it!
      if(!strcmp(parameters.parameters[idx].parameter_name, name)){
        snprintf(parameters.parameters[idx].parameter_value, PARAMETER_LENGTH, "%d", set_value);
        //itoa(set_value, parameters.parameters[idx].parameter_value, 10);
        return TRUE;
      }	
    }
  }
  
  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  snprintf(buffer, PARAMETER_LENGTH, "%d", set_value);
  //  itoa(set_value, buffer, 10);
  result = add_parameter(name, buffer);
    
  return result;
}

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
  )
{
  int idx;
  char *endptr;
  double value;

  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet, using default value");
    return(default_value);
  }
  
  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      /* there is a parameter with the right name.  Now 
	 try to convert it to a double*/
      value = strtod(parameters.parameters[idx].parameter_value, &endptr);
      /*if((value == HUGE_VALF) ||  // AAK removed
	 (value == -HUGE_VALF) || 
	 (endptr == parameters.parameters[idx].parameter_value)) {
	die("Conversion error when trying to convert parameter %s with value %s to an double\n",
	       parameters.parameters[idx].parameter_name,
	    parameters.parameters[idx].parameter_value);*/
      //} else {
	parameters.parameters[idx].used = TRUE;
	return(value);
      // }
    }
  }
  carp(CARP_ERROR, "parameter name: %s doesn't exit, using default value: %.2f", name, default_value);
  return(default_value);
}

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
 )
{
  int idx;
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //convert to string
  sprintf(buffer, "%f", set_value);

  //only check if parameter file has already been parsed
  if(parameter_parsed){
    //check if parameter name already exist?
    for(idx = 0; idx < parameters.num_parameters; idx++){
      //if exist ovewrite it!
      if(!strcmp(parameters.parameters[idx].parameter_name, name)){
        //set to TRUE
        strncpy(parameters.parameters[idx].parameter_value, buffer, PARAMETER_LENGTH);
        return TRUE;
      }	
    }
  }

  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  result = add_parameter(name, buffer);
    
  return result;
}

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, abort.
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter(
  char* name  ///< the name of the parameter looking for -in
  )
{
  int idx;
  char* return_value = NULL;

  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet");
    exit(1);
    return(NULL);
  }
  
  return_value = (char*)mymalloc(sizeof(char) * PARAMETER_LENGTH);

  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      // make sure that there is enough storage allocated in the string
      if((int)strlen(parameters.parameters[idx].parameter_value) 
	 > PARAMETER_LENGTH) {
	die("parameter %s with value %s was too long to copy to string\n",
	    parameters.parameters[idx].parameter_name,
	    parameters.parameters[idx].parameter_value);
      }
      strncpy(return_value,
	      parameters.parameters[idx].parameter_value,
	      PARAMETER_LENGTH);
      parameters.parameters[idx].used = TRUE;
      return(return_value);
    }
  }
  carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
  free(return_value);
  exit(1);
  return(NULL);
}

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is a pointer to the original string
 * Thus, user should no free, good for printing
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter_pointer(
  char* name  ///< the name of the parameter looking for -in
  )
{
  int idx;
    
  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet");
    exit(1);
    return(NULL);
  }
  
  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      return parameters.parameters[idx].parameter_value;
    }
  }
  carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
  exit(1);
  return(NULL);
}

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
 )
{
  int idx;
  BOOLEAN_T result;
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //only check if parameter file has already been parsed
  if(parameter_parsed){
    //check if parameter name already exist?
    for(idx = 0; idx < parameters.num_parameters; idx++){
      //if exist ovewrite it!
      if(!strcmp(parameters.parameters[idx].parameter_name, name)){
        strcpy(parameters.parameters[idx].parameter_value, set_value);
        return TRUE;
      }	
    }
  }

  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  result = add_parameter(name, set_value);
    
  return result;
}

/**
 * Prints the parameters.  If lead_string is not null, preprends it to
 * each line.
 */
void print_parameters(
  char* first_line,  ///< the first line to be printed before the parameter list -in
  char* parameter_filename,  ///< the parameter file name -in
  char* lead_string,  ///< the lead string to be printed before each line -in
  FILE* outstream  ///< the output stream -out
  )
{
  int idx;

  if (lead_string != NULL){
    fprintf(outstream, "%s ",lead_string);
  }
  fprintf(outstream, "%s\n", first_line);

  if (lead_string != NULL){
    fprintf(outstream, "%s ",lead_string);
  }
  fprintf(outstream, "parameter filename: %s\n", parameter_filename);

  if (lead_string != NULL){
    fprintf(outstream, "%s ",lead_string);
  }
  fprintf(outstream, "host=%s date=%s\n", hostname(), date_and_time());
  
  for(idx = 0; idx < parameters.num_parameters; idx++){
    if(lead_string != NULL){
      fprintf(outstream, "%s",lead_string);
    }
    fprintf(outstream, "%s=%s\n",
	    parameters.parameters[idx].parameter_name,
	    parameters.parameters[idx].parameter_value);
  }
}

/**
 * Check to see if any parameters were not used.  Issue a warning.
 */
void check_unused_parameters(void)
{
  int idx;

  for(idx = 0; idx < parameters.num_parameters; idx++){
    if (!parameters.parameters[idx].used) {
      fprintf(stderr, "Warning: Ignoring parameter %s=%s.\n",
	      parameters.parameters[idx].parameter_name,
	      parameters.parameters[idx].parameter_value);
    }
  }
}


/**
 * set the parameters to confirmed, thus blocks any additional changes to the parameters
 */
void parameters_confirmed(){
  parameter_plasticity = FALSE;
}


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
  char* set_value,  ///< the value to be set -in
  BOOLEAN_T required ///< is this a required option -in
  )
{
  int idx;
  
  //check if parameters has been initlialized
  if(!parameter_initialized && !parameter_parsed){
    carp(CARP_ERROR, "must inilialize parameters before copying");
    return FALSE;
  }
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //for required options, there are not in the parameter list, thus must add
  if(required){
    if(add_parameter(name, set_value)){
      return TRUE;
    }
    return FALSE;
  }
  
  //check if parameter name already exist?
  for(idx = 0; idx < parameters.num_parameters; idx++){
    //if exist ovewrite it!
    if(!strcmp(parameters.parameters[idx].parameter_name, name)){
      strcpy(parameters.parameters[idx].parameter_value, set_value);
      return TRUE;
    }	
  }
  
  carp(CARP_ERROR, "incorrect parameter name: %s, from the parameter file", name);
  return FALSE;
} 
