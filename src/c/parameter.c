/******************************************************************************
 * FILE: parameter-file.c
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities.
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
 *\struct stores all the different parameters
 *\brief the array that hold all the parameters
 */
struct parameter_array{
  int num_parameters;   ///< number of the total number of parameters
  struct parameter parameters[NUM_PARAMS]; ///< the paraters
};


/**
 * Global variable
 */

// declare a parameter array to be used
struct parameter_array parameters;
BOOLEAN_T parameter_parsed = FALSE; //have I parsed the parameter file?

/********************************************************************
 *
 * the parameter file is assumed to consist of name/value pairs,
 * separated by an equals sign.
 *
 *******************************************************************/

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

  line = (char*)mycalloc(MAX_LINE_LENGTH, sizeof(char));
  parameters.num_parameters = 0;

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
      parameters.parameters[parameters.num_parameters].used = FALSE;
      strncpy(parameters.parameters[parameters.num_parameters].parameter_name, 
	      line,
	      PARAMETER_LENGTH);
      strncpy(parameters.parameters[parameters.num_parameters].parameter_value,
	      &(line[idx+1]),
	      PARAMETER_LENGTH);
      parameters.num_parameters++;
    }
  }

  // Tell the user what we found.
  if (verbosity > NORMAL_VERBOSE) {
    print_parameters("Using parameters:", parameter_filename, "  ", stderr);
  }

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
  if(!parameter_parsed){
    carp(CARP_WARNING, "parameter file has not been parsed yet, using default value");
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
  return(default_value);
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
  if(!parameter_parsed){
    carp(CARP_WARNING, "parameter file has not been parsed yet, using default value");
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

  return(default_value);
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
  if(!parameter_parsed){
    carp(CARP_WARNING, "parameter file has not been parsed yet, using default value");
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
  return(default_value);
}

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is allocated here and must be freed by the caller.
 * If the value is not found, return NULL.
 * \returns the string value to which matches the parameter name, else returns NULL
 */
char* get_string_parameter(
  char* name  ///< the name of the parameter looking for -in
  )
{
  int idx;
  char* return_value = NULL;

  //check if parameter file has been parsed
  if(!parameter_parsed){
    carp(CARP_WARNING, "parameter file has not been parsed yet, returning NULL");
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
  return(NULL);
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

