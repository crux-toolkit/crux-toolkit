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
#include "hash.h"

/**
 *\struct parameter_hash
 *\brief the hash table that holds all the different parameters
 */
struct parameter_hash{
  int num_parameters;   ///< number of the total number of parameters
  HASH_T* hash; ///< the hash table for parameters
};

/**
 * Global variable
 */

struct parameter_hash parameters_hash_table;
struct parameter_hash* parameters = &parameters_hash_table;
BOOLEAN_T parameter_initialized = FALSE; //have the parameters been initialized?
BOOLEAN_T parameter_parsed = FALSE; //have I parsed the parameter file?
BOOLEAN_T parameter_plasticity = TRUE; //can the parameters be changed?


// parse the parameter file given the filename
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
  
  //allocate the hash table
  parameters->hash = new_hash(NUM_PARAMS);
  
  //set number of parameters to zero
  parameters->num_parameters = 0;

  //set verbosity
  set_int_parameter("verbosity", CARP_ERROR);

  //set parameters
  set_string_parameter("parameter-file", "crux_parameter");
    
  //generate_peptide, create_index parameters  
  set_double_parameter("mass-range", 1);
  set_double_parameter("min-mass", 200);
  set_double_parameter("max-mass", 7200);
  set_int_parameter("max-length", 50);
  set_int_parameter("min-length", 6);
  set_string_parameter("cleavages", "tryptic");
  set_string_parameter("isotopic-mass","average");
  set_string_parameter("redundancy", "redundant");
  set_string_parameter("use-index", "F");
  set_string_parameter("sort", "none");      // mass, length, lexical, none  
  set_boolean_parameter("output-sequence", FALSE);
  set_boolean_parameter("missed-cleavages", FALSE);
  
  //searching peptides
  set_double_parameter("mass-offset", 0);

  //score_peptide_spectrum parameters
  set_double_parameter("beta", 0.075);
  set_double_parameter("max-mz", 4000);
  set_int_parameter("charge", 2);
  set_string_parameter("score-type", "xcorr"); 

  //match_collection parameters
  set_double_parameter("mass-window", 3.0);

  //score_spectrum
  set_string_parameter("prelim-score-type", "sp");
  set_int_parameter("max-rank-preliminary", 500);
  set_int_parameter("max-rank-result", 500);
  set_int_parameter("top-fit-sp", 1000);
  
  //set the top ranking peptides to score for LOGP_*
  set_int_parameter("top-rank-p-value", 1);
  
  //how many peptides to sample for EVD parameter estimation
  set_int_parameter("sample-count", 500);

  //what charge state spectra to run among the ones in ms2 file
  set_string_parameter("spectrum-change", "all");
  set_double_parameter("number-runs", INFINITY);
  
  //match_search
  set_string_parameter("match-output-folder", ".");
  set_string_parameter("output-mode", "binary");
  set_string_parameter("seed", "time");
  set_string_parameter("sqt-output-file", "Prefix of <ms2 input filename>.psm");
  set_double_parameter("spectrum-min-mass", 0.0);
  set_double_parameter("spectrum-max-mass", INFINITY);
  set_int_parameter("top-match", 1);
  set_int_parameter("number-decoy-set", 2);
  
  //match_analysis
  set_string_parameter("algorithm", "percolator");
  set_double_parameter("pi0", 0.9);

  //now we have initialized the parameters
  parameter_initialized = TRUE;
}

/**
 * free heap allocated parameters
 */
void free_parameters(void){
  if(parameter_initialized){
    free_hash(parameters->hash);
  }
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
  // copy the name/value pairs to the right parameter
  ++parameters->num_parameters;
  return add_hash(parameters->hash, 
                  my_copy_string(name), 
                  my_copy_string(set_value));
}

/**
 * copy parameters to parameter hash table
 * there must be a matching name in the parameter hash table
 */
BOOLEAN_T copy_parameter(
  char*     name,  ///< the name of the parameter to add -in
  char* set_value  ///< the value to be added -in                  
  )
{
  
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

  return update_hash_value(parameters->hash, name, set_value);
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
    exit(1);
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
  char *endptr;
  float update_mass;

  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    exit(1);
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
    exit(1);
  }

  while(fgets(line, MAX_LINE_LENGTH, f)==line){
    // if a line begins with "#", it is ignored as a comment 
    // FIXME Also empty lines are ignored, probably amore robust check of new line is required
    if(line[0] != '#' && line[0] != '\n'){
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
	exit(1);
      }
      line[idx] = '\0';
      
      //check if it is amino acid mass update
      if(strlen(line) == 1 && 
         (short int)line[0] >= 'A' && 
         (short int)line[0] <= 'Z'){
        
        update_mass = strtod(&(line[idx+1]), &endptr);
        increase_amino_acid_mass(line[0], update_mass);
      }
      //else, its a parameter value setting
      // copy the name/value pairs to the right parameter
      else if(!copy_parameter(line, &(line[idx+1]))){
        exit(1);
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
 * Each of the following functions searches through the hash table of
 * parameters, looking for one whose name matches the string.  The
 * function returns the corresponding value.
 * \returns TRUE if paramater value is TRUE, else FALSE
 */ 
BOOLEAN_T get_boolean_parameter(
 char*     name  ///< the name of the parameter looking for -in
 )
{
  static char buffer[PARAMETER_LENGTH];
  
  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_ERROR, "parameters has not been set yet");
    exit(1);
  }

  char* value = get_hash_value(parameters->hash, name);

  //can't find parameter
  if(value == NULL){
    carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
  
  // make sure that there is enough storage allocated in the string
  if((int)strlen(value) 
     > PARAMETER_LENGTH) {
    die("parameter %s with value %s was too long to copy to string\n",
        name,
        value);
  }
  strncpy(buffer,
          value,
          PARAMETER_LENGTH);
  if (strcmp(buffer, "TRUE") == 0) {
    return(TRUE);
  } 
  else if (strcmp(buffer, "FALSE") == 0) {
    return(FALSE);
  } 
  else {
    die("Invalid Boolean parameter %s.\n", buffer);
  }
  
  carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
  exit(1);
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
  BOOLEAN_T result;
    
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }

  //only check if parameter file has already been parsed
  if(parameter_parsed || parameter_initialized){
    if(set_value){
      return update_hash_value(parameters->hash, name, "TRUE");
    }
    else{
      return update_hash_value(parameters->hash, name, "FALSE");
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
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the int value of the parameter
 */
int get_int_parameter(
  char* name  ///< the name of the parameter looking for -in
  )
{
  char *endptr;
  long int value;

  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_ERROR, "parameters has not been set yet");
    exit(1);
  }

  char* int_value = get_hash_value(parameters->hash, name);

  //can't find parameter
  if(int_value == NULL){
    carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
  
  /* there is a parameter with the right name.  Now 
     try to convert it to a base 10 integer*/
  value = strtol(int_value, &endptr, 10);
  if ((value == LONG_MIN) || 
      (value == LONG_MAX) || 
      (endptr == int_value)) {
    die("Conversion error when trying to convert parameter %s with value %s to an int\n",
        name, 
        int_value);
  } 
  
  return((int)value);
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
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //only check if parameter file has already been parsed
  if(parameter_parsed  || parameter_initialized){
    snprintf(buffer, PARAMETER_LENGTH, "%d", set_value);
    return update_hash_value(parameters->hash, name, buffer);
  }
  
  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  snprintf(buffer, PARAMETER_LENGTH, "%d", set_value);  
  result = add_parameter(name, buffer);
  
  return result;
}

/**
 * Searches through the list of parameters, looking for one whose
 * name matches the string.  This function returns the parameter value if the
 * parameter is in the parameter hash table.  This
 * function exits if there is a conversion error. 
 *\returns the double value of the parameter
 */
double get_double_parameter(
  char* name   ///< the name of the parameter looking for -in
  )
{
  char *endptr;
  double value;
  
  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_ERROR, "parameters has not been set yet");
    exit(1);
  }

  char* double_value = get_hash_value(parameters->hash, name);
 
  //can't find parameter
  if(double_value == NULL){
    carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
 
  /* there is a parameter with the right name.  Now 
     try to convert it to a double*/
  value = strtod(double_value, &endptr);
  /*if((value == HUGE_VALF) ||  // AAK removed
    (value == -HUGE_VALF) || 
    (endptr == double_value)) {
    die("Conversion error when trying to convert parameter %s with value %s to an double\n",
    name,
    double_value);*/
  //} else {  
  return(value);
  // }
  
  carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
  exit(1);
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
  if(parameter_parsed || parameter_initialized){
    return update_hash_value(parameters->hash, name, buffer);    
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
  
  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet");
    exit(1);
    return(NULL);
  }
  
  char* string_value = get_hash_value(parameters->hash, name);
  
  //can't find parameter
  if(string_value == NULL){
    carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
  
  return my_copy_string(string_value);
}

/**
 * Searches through the list of parameters, looking for one whose
 * parameter_name matches the string. 
 * The return value is a pointer to the original string
 * Thus, user should not free, good for printing
 * \returns the string value to which matches the parameter name, else aborts
 */
char* get_string_parameter_pointer(
  char* name  ///< the name of the parameter looking for -in
  )
{
  //check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet");
    exit(1);
    return(NULL);
  }
  
  char* string_value = get_hash_value(parameters->hash, name);

  //can't find parameter
  if(string_value == NULL){
    carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
  else{
    return string_value;
  }
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
  BOOLEAN_T result;
  
  //check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //only check if parameter file has already been parsed
  if(parameter_parsed  || parameter_initialized){
    return update_hash_value(parameters->hash, name, set_value);
  }
  
  //if it doesn't already exist(wasn't in the parameter file), add to parameter list
  result = add_parameter(name, set_value);
    
  return result;
}

/**
 * Prints the parameters.  If lead_string is not null, preprends it to
 * each line.
 */
/*
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
*/
/**
 * set the parameters to confirmed, thus blocks any additional changes to the parameters
 */
void parameters_confirmed(){
  parameter_plasticity = FALSE;
}


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
  char* set_value,  ///< the value to be set -in
  BOOLEAN_T required ///< is this a required option -in
  )
{
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

  //if exist ovewrite it!
  return update_hash_value(parameters->hash, name, set_value);  
} 

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

