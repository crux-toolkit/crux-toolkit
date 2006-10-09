/******************************************************************************
 * FILE: parameter-file.c
 * AUTHOR: Tobias Mann and Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities.
 *****************************************************************************/
#include "utils.h"
#include "parameter-file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

// The parameters are stored as a global variable.
struct parameter{
  BOOLEAN_T used;
  char parameter_name[PARAMETER_LENGTH];
  char parameter_value[PARAMETER_LENGTH];
};

struct parameter_array{
  int num_parameters;
  struct parameter parameters[NUM_PARAMS];
};

struct parameter_array parameters;

/********************************************************************
 *
 * the parameter file is assumed to consist of name/value pairs,
 * separated by an equals sign.
 *
 *******************************************************************/
void parse_parameter_file(char* parameter_filename){
  FILE *f;
  char *line;
  int idx;

  line = (char*) calloc(MAX_LINE_LENGTH, sizeof(char));
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
}

/*******************************************************************
 * 
 * Each of the following functions searches through the list of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value, or the default value if
 * the parameter is not found.
 *
 *******************************************************************/
BOOLEAN_T get_boolean_parameter
(char*     name, 
 BOOLEAN_T default_value)
{
  int idx;
  static char buffer[PARAMETER_LENGTH];

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

/*******************************************************************
 * 
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string.  This function returns 1 if the
 * parameter is in the parameter_array, and 0 otherwise.  This
 * function exits if there is a conversion error. If all goes well,
 * then the value of the parameter is written to the passed pointer i.
 * If there is an error, then the value of *i is not changed.
 *
 *******************************************************************/
int get_int_parameter
(char* name, 
 int default_value)
{
  int idx;
  char *endptr;
  long int value;

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

/*******************************************************************
 * 
 * Same as above, but for a double.
 *
 *******************************************************************/
double get_double_parameter
(char* name,
 double default_value)
{
  int idx;
  char *endptr;
  double value;

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

/*******************************************************************
 * 
 * Similar to above, but for a string.  The return value is allocated
 * here and must be freed by the caller.  If the value is not found,
 * return NULL.
 *
 *******************************************************************/
char* get_string_parameter
(char* name)
{
  int idx;
  char* return_value = NULL;

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

/**********************************************************************
 *
 * Prints the parameters.  If lead_string is not null, preprends it to
 * each line.
 *
 *********************************************************************/
void print_parameters
(char* first_line,
 char* parameter_filename,
 char* lead_string,
 FILE* outstream)
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

/**********************************************************************
 *
 * Check to see if any parameters were not used.  Issue a warning.
 *
 *********************************************************************/
void check_unused_parameters ()
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

