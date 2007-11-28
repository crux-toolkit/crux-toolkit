/****************************************************************************
 * FILE: parameter.c
 * AUTHOR: written by Tobias Mann, CRUXified by Chris Park
 * CREATE DATE: 2006 Oct 09
 * DESCRIPTION: General parameter handling utilities. MUST declare ALL optional command parameters here inside initalialize_parameters
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "spectrum.h"
#include "peak.h"
#include "carp.h"
#include "mass.h"
#include "scorer.h"
#include "utils.h"
#include "parameter.h"
#include "parse_arguments.h"
#include "hash.h"
#define MAX_SET_PARAMS 256
//TODO:  why are #includes here and not in .h?
//       change all temp_add to use add_or_update and no strcpy
//       in all temp_set, change result=add_... to result= result && add_...

/**
 *\struct parameter_hash
 *\brief the hash table that holds all the different parameters
 */
struct parameter_hash{
  int num_parameters;   ///< number of the total number of parameters
  HASH_T* hash; ///< the hash table for parameters
};

/**
 * Global variables
 */
static char* parameter_type_strings[NUMBER_PARAMETER_TYPES] = { "INT_ARG", "DOUBLE_ARG", "STRING_ARG", "MASS_TYPE_T", "PEPTIDE_TYPE_T"};

//one hash for parameter values, one for usage statements, one for types
struct parameter_hash  parameters_hash_table;
struct parameter_hash* parameters = &parameters_hash_table;
struct parameter_hash  usage_hash_table;
struct parameter_hash* usages = &usage_hash_table;
struct parameter_hash  type_hash_table;
struct parameter_hash* types = & type_hash_table;
/* MINMAX
struct parameter_hash  min_values_hash_table;
struct parameter_hash* min_values = & type_hash_table;
struct parameter_hash  max_values_hash_table;
struct parameter_hash* max_values = & type_hash_table;
*/
char* parameters_set[MAX_SET_PARAMS]; // list of option names set on cmd line or in param file
int num_params_set = 0;

BOOLEAN_T parameter_initialized = FALSE; // have the parameters been initialized?
BOOLEAN_T usage_initialized = FALSE; // have the usages been initialized?
BOOLEAN_T type_initialized = FALSE; // have the types been initialized?

//remove this?
BOOLEAN_T parameter_parsed = FALSE; // have I parsed the parameter file?
BOOLEAN_T parameter_plasticity = TRUE; // can the parameters be changed?

/************************************
 * Private function declarations
 ************************************ 
 */

// parse the parameter file given the filename
// called by parse command line
void parse_parameter_file(
  char* parameter_filename ///< the parameter file to be parsed -in
  );

/**
 * Examine each option in list to determine if the values
 * are within the proper range and of the correct type
 * Requires (or at least only makes sense after) 
 * parse_cmd_line_into_params_hash() has been run.
 */
BOOLEAN_T check_option_type_and_bounds(char* name);

/**
 *
 */
BOOLEAN_T string_to_param_type(char*, PARAMETER_TYPE_T* );

/**
 * temporary replacement for function, return name once all exe's are fixed
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T temp_set_boolean_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T set_value,  ///< the value to be set -in
 char* usage
 );

BOOLEAN_T temp_set_int_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the value to be set -in
 int max_value,  ///< the value to be set -in
 char* usage
 );

BOOLEAN_T temp_set_double_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 char* usage
  );

BOOLEAN_T temp_set_string_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 char* set_value,  ///< the value to be set -in
 char* usage
  );

BOOLEAN_T temp_set_mass_type_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 char* usage
  );

BOOLEAN_T temp_set_peptide_type_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 PEPTIDE_TYPE_T set_value,  ///< the value to be set -in
 char* usage
  );

BOOLEAN_T select_cmd_line(  
  char** option_names, ///< list of options to be allowed for main -in
  int    num_options,  ///< number of optons in that list -in
  int (*parse_argument_set)(char*, char*, void*, enum argument_type) ///< function point to choose arguments or options 
  );

  //double get_option_min(char*);
  //double get_option_max(char*);
//are there string parameters that have a limited number
// of accepted strings that are not also a defined type?
//BOOLEAN_T get_option_accepted_list(char* name, char** put_list_here, int& num);

/************************************
 * Function definitions
 ************************************
 */


/**
 * initialize parameters
 * ONLY add optional parameters here!!!
 * MUST declare ALL optional parameters in array to be used
 * Every option and its default value for every executable 
 * must be declared here
 */
void initialize_parameters(void){
  carp(CARP_DETAILED_DEBUG, "Initializing parameters in parameter.c");

  // check if parameters been initialized
  if(parameter_initialized){
    carp(CARP_ERROR, "parameters have already been initialized");
    return;
  }
  
  // allocate the hash tables
  parameters->hash = new_hash(NUM_PARAMS);
  usages->hash = new_hash(NUM_PARAMS);
  types->hash = new_hash(NUM_PARAMS);
  /* MINMAX
  min_values->hash = new_hash(NUM_PARAMS);
  max_values->hash = new_hash(NUM_PARAMS);
  */  

  // set number of parameters to zero
  parameters->num_parameters = 0;
  usages->num_parameters = 0;
  types->num_parameters = 0;
  /* MINMAX
  min_values->num_parameters = 0;
  max_values->num_parameters = 0;
  */

  // set verbosity
  temp_set_int_parameter("verbosity", CARP_ERROR, CARP_FATAL, CARP_MAX,
	"Set level of output to stderr (0-100).  Default 50.");

  // set parameter file name (no default)
  //set_string_parameter("parameter-file", "crux.params");
  temp_set_string_parameter("parameter-file", NULL, 
	"Set additional options with values in the given file.");
    
  // generate_peptide arguments
  temp_set_string_parameter("protein input", NULL, 
  "File containing protein sequences either in fasta format or binary index.");
  // create_index arguments
  temp_set_string_parameter("protein fasta file", NULL,
		    "File containing protein sequences in fasta format.");
  temp_set_string_parameter("index name", NULL,
		    "Name to give the new directory containing index files.");

  // generate_peptide, create_index parameters  
  temp_set_double_parameter("min-mass", 200, 0, 7200,
	"The minimum mass of peptides to consider. Default 200.");
  temp_set_double_parameter("max-mass", 7200, 1, BILLION, 
	"The maximum mass of peptides to consider. Default 7200.");
  temp_set_int_parameter("min-length", 6, 1, MAX_PEPTIDE_LENGTH,
	"The minimum length of peptides to consider. Default 6.");
  temp_set_int_parameter("max-length", 50, 1, MAX_PEPTIDE_LENGTH,
	"The maximum length of peptides to consider. Default 50.");
  temp_set_boolean_parameter("missed-cleavages", FALSE, 
	"Include peptides with missed cleavage sites. Default FALSE.");
  temp_set_peptide_type_parameter("cleavages", TRYPTIC, 
	"The type of cleavage sites to consider (tryptic, partial, all)");
  temp_set_mass_type_parameter("isotopic-mass", AVERAGE, 
	"Which isotopes to use in calcuating mass (average or mono). " \
	"Default average");
  
  // more generate_peptide parameters
  temp_set_boolean_parameter("output-sequence", FALSE, "usage");
  temp_set_boolean_parameter("output-trypticity", FALSE, "usage");
  temp_set_string_parameter("use-index", "F", "usage");
  temp_set_string_parameter("sort", "none", "usage");//mass,length,lexical,none  
  temp_set_boolean_parameter("unique-peptides", FALSE, "usage");

  // more create_index parameters
  //temp_set_double_parameter("mass-range", 10, "usage");

    // searching peptides
  temp_set_double_parameter("mass-offset", 0.0, 0, 0, "usage");

  // score_peptide_spectrum parameters
  temp_set_double_parameter("beta", 0.075, 0, 1, "usage");
  temp_set_double_parameter("max-mz", 4000, 0, BILLION, "usage");
  temp_set_int_parameter("charge", 2, 1, 4, "usage");
  temp_set_string_parameter("score-type", "xcorr", "usage"); 

  // match_collection parameters
  temp_set_double_parameter("mass-window", 3.0, 0, 100, "usage");

  // create_psm_files
  temp_set_int_parameter("starting-sentence-idx", 0, 0, 0, "usage");
  //set_string_parameter("model-type", "single");
  temp_set_string_parameter("model-type", "single", "usage");

  // score_spectrum
  temp_set_string_parameter("prelim-score-type", "sp", "usage");
  temp_set_int_parameter("max-rank-preliminary", 500, 1, BILLION, "usage");
  temp_set_int_parameter("max-rank-result", 500, 1, BILLION, "usage");
  temp_set_int_parameter("top-fit-sp", 1000, 1, BILLION, "usage");
  temp_set_int_parameter("number-top-scores-to-fit", -1, -10, BILLION, "usage");
  temp_set_int_parameter("number-peptides-to-subset", 0, 0, 0, "usage");
  temp_set_double_parameter("fraction-top-scores-to-fit", -1.0, -10, 10, "usage");
  temp_set_int_parameter("skip-first-score", 0, 0, 1, "usage");
  
  // set the top ranking peptides to score for LOGP_*
  temp_set_int_parameter("top-rank-p-value", 1, 1, BILLION, "usage");
  
  // how many peptides to sample for EVD parameter estimation
  temp_set_int_parameter("sample-count", 500, 0, BILLION, "usage");

  // what charge state spectra to run among the ones in ms2 file
  temp_set_string_parameter("spectrum-charge", "all", "usage");
  temp_set_double_parameter("number-runs", BILLION, 1, BILLION, "usage");
  
  // match_search
  temp_set_string_parameter("match-output-folder", ".", "usage");
  temp_set_string_parameter("output-mode", "binary", "usage"); //binary, sqt, all
  temp_set_string_parameter("seed", "time", "usage");
  temp_set_string_parameter("sqt-output-file", "target.psm", "usage");
  temp_set_string_parameter("decoy-sqt-output-file", "decoy.psm", "usage");
  temp_set_double_parameter("spectrum-min-mass", 0.0, 0, BILLION, "usage");
  temp_set_double_parameter("spectrum-max-mass", BILLION, 1, BILLION, "usage");
  temp_set_int_parameter("top-match", 1, 1, 111, "usage");
  temp_set_int_parameter("number-decoy-set", 2, 0, 10, "usage");
  
  // match_analysis
  temp_set_string_parameter("algorithm", "percolator", "usage");
  temp_set_string_parameter("feature-file", "match_analysis.features", "usage");
  temp_set_double_parameter("pi0", 0.9, 0, 1, "usage");
  temp_set_string_parameter("percolator-intraset-features", "F", "usage"); // for false

  //static mods
  /*  The are commented out only because I am too lazy to give min and max
  temp_set_double_parameter("A", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("B", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("C", 57.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("D", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("E", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("F", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("G", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("H", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("I", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("J", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("K", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("L", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("M", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("N", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("O", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("P", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("Q", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("R", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("S", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("T", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("U", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("V", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("W", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("X", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("Y", 0.0, "NOT FOR COMMAND LINE");
  temp_set_double_parameter("Z", 0.0, "NOT FOR COMMAND LINE");
  */

  // now we have initialized the parameters
  parameter_initialized = TRUE;
  usage_initialized = TRUE;
  type_initialized = TRUE;


}


/*
 * Main calls this to determine which required arguments
 * must be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_req()
 */
BOOLEAN_T select_cmd_line_arguments(  //remove options from name
  char** option_names,
  int    num_options //, int (*parse_argument_set)(char, char, void, enum argument_type) 
  ){
  select_cmd_line( option_names, num_options, 
		   &parse_arguments_set_req);
  return TRUE;
}

/*
 * Main calls this to determine which options
 * can be used.  Must be called after initialize_parameters
 * This is the interface for parse_argument_set_opt()
 */
BOOLEAN_T select_cmd_line_options(  //remove options from name
  char** option_names,
  int    num_options //, int (*parse_argument_set)(char, char, void, enum argument_type) 
  ){
  select_cmd_line( option_names, num_options, 
		   &parse_arguments_set_opt);
  return TRUE;
}
/*
 * Private function for doing the work of select_cmd_line_options
 * and select_cmd_line_arguments which is all the same except for
 * the last function call which is now set with a function pointer
 * 
 */
BOOLEAN_T select_cmd_line(  //remove options from name
  char** option_names,
  int    num_options, 
  int (*parse_arguments_set_ptr)(char*, char*, void*, enum argument_type) 
  ){

  carp(CARP_DETAILED_DEBUG, "Selecting options");
  BOOLEAN_T success = TRUE;

  if( (num_options < 1) || (option_names == NULL) ){
    success = FALSE; //?
    return success;
  }

  //for each name in list
  int i;
  for( i=0; i< num_options; i++){
    carp(CARP_DETAILED_DEBUG, "Option is: %s", option_names[i]);
    //get value, usage, types
    void* value_ptr = get_hash_value(parameters->hash, option_names[i]);
    void* usage_ptr = get_hash_value(usages->hash, option_names[i]);
    void* type_ptr =  get_hash_value(types->hash, option_names[i]);
    if( strcmp(type_ptr, "PEPTIDE_TYPE_T") == 0 ||
	strcmp(type_ptr, "MASS_TYPE_T") == 0){
      type_ptr = "STRING_ARG";
    }
    carp(CARP_DETAILED_DEBUG, "Found value: %s, usage: %s, type: %s", (char*)value_ptr, (char*)usage_ptr, (char*)type_ptr);


    // check that it is in the params hash
    // since default value can be null, don't check value
    if( value_ptr == NULL || usage_ptr == NULL || type_ptr == NULL ){
      carp(CARP_FATAL, 
	   "Cannot select parameter '%s'. Value, usage or type not found.\nFound value: %s, usage: %s, type: %s", 
	   option_names[i],
	   value_ptr,
	   usage_ptr,
	   type_ptr);
      
      exit(1);  // or  set success to F?
    }
    // add the option via parse_arguments    
    success = parse_arguments_set_ptr(option_names[i],
				      usage_ptr,
				      value_ptr, 
				      string_to_argument_type(type_ptr)); 
  }

  carp(CARP_DETAILED_DEBUG, "Did setting the arguments work? %i", success);
  return success;
}
/**
 * helper used below.  look for param file name, die if error
 * return null if not found
 */
BOOLEAN_T find_param_filename(int argc, 
			      char** argv, 
			      char* filename_buffer, 
			      int buffer_size){
  BOOLEAN_T success = TRUE;
  int i;
  int param_file_index = -1;
  for( i=0; i< argc; i++){
    if( strcmp(argv[i], "--parameter-file") == 0){
      param_file_index = i+1;
      break;
    }
  }
  // check for error
  if( param_file_index >= argc ){
    carp(CARP_FATAL, "Option '--parameter-file' requires argument");
    exit(1);
  }
  //return the filename
  else if( param_file_index > 0 ){  
    char* param_filename = argv[param_file_index];
    carp(CARP_DETAILED_DEBUG, "Parameter file name is %s", param_filename);

    if( strlen(param_filename) < (unsigned)buffer_size ){
      strcpy(filename_buffer, param_filename);
      success = TRUE;
    }
    else{
      carp(CARP_FATAL, "Parameter filename is too long");
      exit(1);
    }
  }
  else{ //parameter_file_index < 0, i.e. no paramter file option
    success = FALSE;
  }

  return success;
}
/**
 * Take the command line string from main, find the parameter fil
 * option (if present), parse it's values into the hash, and parse
 * the command line options and arguments into the hash
 * main then retrieves the values through get_value
 */
BOOLEAN_T parse_cmd_line_into_params_hash(int argc, char** argv){
  carp(CARP_DETAILED_DEBUG, "Parameter.c is parsing the command line");
  BOOLEAN_T success = TRUE;
  int i;
  //first look for parameter-file option and parse those values before
  //    command line values

  // search argv for parameter file
  // TODO find appropriate const for 256
  char param_filename[256];
  if(find_param_filename(argc, argv, param_filename, 256)){
    parse_parameter_file(param_filename);  //this checks types and bounds
  }
  else{ 
    carp(CARP_INFO, "No parameter file specified.  Using defaults and command line values");
  }

  // now parse the command line and put those values in hash
  //  overwriting file parameters

  success = parse_arguments_into_hash(argc, argv, parameters->hash, 0); 
  if( success ){
    // check each option value
    for(i=1; i<argc; i++){
      char* word = argv[i];
      if( word[0] == '-' ){   //if word starts with --
	word = word + 2;      //ignore the --
	check_option_type_and_bounds(word);
      }//else skip this word
    }

  }
  else{  // parse_arguments encountered an error
    char* error_message = NULL;
    char* usage = parse_arguments_get_usage("create_index");
    int error_code = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", error_code);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
    exit(1);
  }
  
  return success;
}

/*
 * this will mimic initialize_parameters, having knowledge of all
 * option names, types and allowed bounds
 * or maybe it will be generalized and there will be yet another
 * hash with min and max values.
 */
BOOLEAN_T check_option_type_and_bounds(char* name){

  BOOLEAN_T success = TRUE;
  char* die_str;
  char* type_str = get_hash_value(types->hash, name);
  char* value_str = get_hash_value(parameters->hash, name);
  /* MINMAX
  char* min_str = get_hash_value(min_values->hash, name);
  char* max_str = get_hash_value(max_values->hash, name);
  int min_int;
  int max_int;
  double value_double;
  double min_double;
  double max_double;
  */

  MASS_TYPE_T mass_type;
  PEPTIDE_TYPE_T pep_type;

  PARAMETER_TYPE_T param_type;
  string_to_param_type( type_str, &param_type ); 

  carp(CARP_DETAILED_DEBUG, 
       "Checking option '%s' of type '%s' for type and bounds", 
       name, type_str);

  switch( param_type ){
  case INT_P:
    carp(CARP_DETAILED_DEBUG, "found int opt with value %i\n", 
	 atoi(value_str));
    //atoi, check min/max
    /* MINMAX
    value_int = atoi(value_str);
    min_int = atoi(min_str);
    max_int = atoi(max_str);
    printf("found int opt with value %d, min %d, max %d\n", value_int, min_int, max_int);
    */
    break;
  case DOUBLE_P:
     carp(CARP_DETAILED_DEBUG, "found double opt with value %f\n", 
	  atof(value_str));
    //atof, check min/max
    /* MINMAX
    value_double = atof(value_str);
    min_double = atof(min_str);
    max_double = atof(max_str);
    printf("found double opt with value %f, min %f, max %f\n", value_double, min_double, max_double);
    */
    break;
  case STRING_P:
    carp(CARP_DETAILED_DEBUG, "found string opt with value %s\n", value_str);
    //check list of legal values?
    break;
  case MASS_TYPE_P:
    carp(CARP_DETAILED_DEBUG, "found mass_type opt with value %s\n", 
	 value_str);
    if( ! string_to_mass_type( value_str, &mass_type )){
      success = FALSE;
      die_str = "Illegal mass-type.  Must be 'mono' or 'average'";
    }
    break;
  case PEPTIDE_TYPE_P:
      carp(CARP_DETAILED_DEBUG, "found peptide_type param, value '%s'\n", 
	   value_str);
    if( ! string_to_peptide_type( value_str, &pep_type )){
      success = FALSE;
      die_str = "Illegal peptide cleavages.  Must be...something";
    }
    break;
  default:
    carp(CARP_FATAL, "Your param type wasn't found");
    exit(1);
  }

  if( ! success ){
    carp(CARP_FATAL, die_str);
    exit(1);
  }
  return success;
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
  
  // check if parameters has been initlialized
  if(!parameter_initialized){
    carp(CARP_ERROR, "must inilialize parameters before copying");
    return FALSE;
  }
    
  // check if parameters can be changed
  if(!parameter_plasticity){
   carp(CARP_ERROR, "can't change parameters once they are confirmed");
   return FALSE;
  }
  //how does memory work for update? allocate new?
  return update_hash_value(parameters->hash, name, set_value);
}

/**
 * This method should be called only after parsed command line
 * first, parse paramter file
 * Next, update the parameter files with command line options
 * command line arguments have higher precedence 
 * parse the parameter file given the filename
 */
void parse_update_parameters(
  char* parameter_file ///< the parameter file to be parsed -in
  )
{
  // initialize
  initialize_parameters();
  
  // if no parameter file name has been specified in command line,
  // use default parameter filename
  if(parameter_file != NULL){
    // parse parameter file
    parse_parameter_file(parameter_file);
  }
  else{
    // parse parameter file
    //parse_parameter_file(get_string_parameter_pointer("parameter-file"));
  }
  
  // update the parameters if any comman line arguments exist
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
  FILE *file;
  char *line;
  int idx;
  //  char *endptr;
  //  float update_mass;

  carp(CARP_DETAILED_DEBUG, "Parsing parameter file '%s'",parameter_filename);

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_FATAL, "Can't change parameters once they are confirmed");
    exit(1);
  }

  // check if parameter file exists, if not die 
  if(access(parameter_filename, F_OK)){
    carp(CARP_FATAL, "Could not open parameter file.");
    exit(1);
  }

  file = fopen(parameter_filename, "r");
  if(file == NULL){
    //change to stderr
    carp(CARP_FATAL, "Couldn't open parameter file '%s'", parameter_filename);
    exit(1);
  }

  line = (char*)mycalloc(MAX_LINE_LENGTH, sizeof(char));

  while(fgets(line, MAX_LINE_LENGTH, file)==line){
    //    printf("read line '%s'.\n", line);

    idx = 0;
    
    // Change the newline to a '\0' ignoring trailing whitespace
    for(idx = MAX_LINE_LENGTH - 1; idx >= 0; idx--){
      if(line[idx] == '\n' || line[idx] == '\r' || 
	 line[idx] == '\f' || line[idx] == ' ' || line[idx] == '\t')
	line[idx] = '\0';
    }
    /* why does this segfault?  only with break, not without
    if(line[0] == '#' || line[0] == '\0'){
      printf("comment or blank line");
      break;
    }
    */
    // empty lines and those beginning with '#' are ignored
    if(line[0] != '#' && line[0] != '\0'){

      /* find the '=' in the line.  Exit with error if the line 
	 has no equals sign. */
      while(idx < (int)strlen(line) && line[idx] != '='){
	idx++;
      }
      if(idx == 0 || idx >= (int)(strlen(line)-1)){
	//these should be FATALs
	carp(CARP_ERROR, "Lines in a parameter file must have the form:\n");
	carp(CARP_ERROR, "\n\tname=value\n\n");
	carp(CARP_ERROR, 
	     "In file %s, the line\n%s\ndoes not have this format\n",
	     parameter_filename, line);
	exit(1);
      }

      line[idx] = '\0';
      char* option_name = line;
      char* option_value = &(line[idx+1]);
      carp(CARP_DETAILED_DEBUG, "Found option '%s' and value '%s'", option_name, option_value);

      if(! update_hash_value(parameters->hash, option_name, option_value) ){
	carp(CARP_ERROR, "Unexpected parameter file option '%s'", option_name);
	exit(1);
      }

      //      parameters_set[num_params_set++] = my_copy_string(option_name);
      check_option_type_and_bounds(option_name);

      /*
      // check if it is amino acid mass update
      if(strlen(line) == 1 && 
         (short int)line[0] >= 'A' && 
         (short int)line[0] <= 'Z'){
        
        update_mass = strtod(&(line[idx+1]), &endptr);
        increase_amino_acid_mass(line[0], update_mass);
      }
      // else, its a parameter value setting
      // copy the name/value pairs to the right parameter
      else if(!copy_parameter(line, &(line[idx+1]))){
        exit(1);
      }
      */
    }
  }

  fclose(file);
  myfree(line);

  // now we have parsed the parameter file
  parameter_parsed = TRUE;
}

/**************************************************
 *   GETTERS (public)
 **************************************************
 */

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
  
  // check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_ERROR, "parameters has not been set yet");
    exit(1);
  }

  char* value = get_hash_value(parameters->hash, name);

  // can't find parameter
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

  if ((strcmp(buffer, "TRUE") == 0) || (strcmp(buffer, "T") == 0)){
    return(TRUE);
  } 
  else if ((strcmp(buffer, "FALSE") == 0) || (strcmp(buffer, "F") == 0)){
    return(FALSE);
  } 
  else {
    die("Invalid Boolean parameter %s.\n", buffer);
  }
  
  carp(CARP_FATAL, "parameter name: %s, doesn't exit", name);
  exit(1);
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

  // check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_ERROR, "parameters has not been set yet");
    exit(1);
  }

  char* int_value = get_hash_value(parameters->hash, name);

  // can't find parameter
  if(int_value == NULL){
    carp(CARP_FATAL, "parameter name: %s, doesn't exit", name);
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
  
  // check if parameter file has been parsed
  if(!parameter_initialized){
    carp(CARP_FATAL, "parameters have not been set yet");
    exit(1);
  }

  char* double_value = get_hash_value(parameters->hash, name);
 
  // can't find parameter
  if(double_value == NULL){
    carp(CARP_FATAL, "parameter name '%s', doesn't exit", name);
    exit(1);
  }
 
  /* there is a parameter with the right name.  Now 
     try to convert it to a double*/
  value = strtod(double_value, &endptr);
  /*if((value == HUGE_VALF) ||  // AAK removed //BF: why?
    (value == -HUGE_VALF) || 
    (endptr == double_value)) {
    die("Conversion error when trying to convert parameter %s with value %s to an double\n",
    name,
    double_value);*/
  // } else {  
  return(value);
  // }
  
  carp(CARP_ERROR, "parameter name: %s, doesn't exit", name);
  exit(1);
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
  
  // check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_WARNING, "parameters has not been set yet");
    exit(1);
    return(NULL);
  }
  
  char* string_value = get_hash_value(parameters->hash, name);
  
  // can't find parameter
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
  // check if parameter file has been parsed
  if(!parameter_parsed && !parameter_initialized){
    carp(CARP_FATAL, "parameters has not been set yet");
    exit(1);
  }
  
  char* string_value = get_hash_value(parameters->hash, name);

  // can't find parameter
  if(string_value == NULL){
    carp(CARP_FATAL, "parameter name: %s, doesn't exit", name);
    exit(1);
  }
  else{
    return string_value;
  }
}

PEPTIDE_TYPE_T get_peptide_type_parameter(
  char* name
    ){

  char* param = get_string_parameter_pointer(name);
  /*
  int peptide_type = convert_enum_type_str(
      param, 0, peptide_type_strings, NUMBER_PEPTIDE_TYPES);
  */
  PEPTIDE_TYPE_T peptide_type;
  int success = string_to_peptide_type(param, &peptide_type);
  //we should have already checked the type, but just in case
  if( success < 0 ){
    carp(CARP_FATAL, "Peptide_type parameter %s has the value %s which is not of the correct type\n", name, param);
    exit(1);
  }
  return peptide_type;
}

MASS_TYPE_T get_mass_type_parameter(
   char* name
   ){
  char* param_value_str = get_hash_value(parameters->hash, name);
  MASS_TYPE_T param_value;
  BOOLEAN_T success = string_to_mass_type(param_value_str, &param_value);

  if( ! success ){
    carp(CARP_FATAL, "Mass_type parameter %s has the value %s which is not of the correct type", name, param_value_str);
    exit(1);
  }
  return param_value;
}

/**************************************************
 *   SETTERS (private)
 **************************************************
 */

//these could all use add instead of add_or_update
BOOLEAN_T temp_set_boolean_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 BOOLEAN_T set_value,  ///< the value to be set -in
 char* usage
  )
{
  BOOLEAN_T result;
    
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }

  char* bool_str;
  if(set_value){
    bool_str = "TRUE";
  }
  else{
    bool_str = "FALSE";
  }
  result = add_or_update_hash(parameters->hash, name, bool_str);
  result = add_or_update_hash(usages->hash, name, usage);
  result = add_or_update_hash(types->hash, name, "BOOL_ARG");

  return result;
}

//temporary, replace name with set_int_parameter
BOOLEAN_T temp_set_int_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 int set_value,  ///< the value to be set -in
 int min_value,  ///< the value to be set -in
 int max_value,  ///< the value to be set -in
 char* usage
  )
{
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  //stringify default, min, and max values and set
  snprintf(buffer, PARAMETER_LENGTH, "%i", set_value);
  result = add_or_update_hash(parameters->hash, name, buffer);

   carp(CARP_DETAILED_DEBUG, "not setting min %i or max %i\n", 
	min_value, max_value);
  /* MINMAX
  snprintf(buffer, PARAMETER_LENGTH, "%i", min_value);
  result = add_or_update_hash(min_values->hash, name, buffer);

  snprintf(buffer, PARAMETER_LENGTH, "%i", max_value);
  result = add_or_update_hash(max_values->hash, name, buffer);
  */

  result = add_or_update_hash(usages->hash, name, usage);
  result = add_or_update_hash(types->hash, name, "INT_ARG");
  
  return result;
}


//change name when all exe's are fixed
BOOLEAN_T temp_set_double_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 char* usage
  )
{
  BOOLEAN_T result;
  char buffer[PARAMETER_LENGTH];
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // convert to string
  snprintf(buffer, PARAMETER_LENGTH, "%f", set_value);
  result = add_or_update_hash(parameters->hash, name, buffer);    

   carp(CARP_DETAILED_DEBUG, "not setting min %f max %f\n", 
	min_value, max_value);
  /* MINMAX
  snprintf(buffer, PARAMETER_LENGTH, "%f", min_value);
  result = add_or_update_hash(min_values->hash, name, buffer);    
  printf("\tsetting double min %s", buffer);

  snprintf(buffer_max, PARAMETER_LENGTH, "%f", max_value);
  result = add_or_update_hash(max_values->hash, name, buffer_max);    
  printf("\tsetting double max %s\n", buffer_max);

  temp = get_hash_value(max_values->hash, name);
  printf("Got back out max of %s\n", temp);
  */

  result = add_or_update_hash(usages->hash, name, usage);    
  result = add_or_update_hash(types->hash, name, "DOUBLE_ARG");    

  //  printf("Got back out max of %s\n", temp);  //MINMAX
  return result;
}
/**
 * temporary replacement for function, return name once all exe's are fixed
 * \returns TRUE if paramater value is set, else FALSE
 */ 
BOOLEAN_T temp_set_string_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 char* set_value,  ///< the value to be set -in
 char* usage
  )
{
  BOOLEAN_T result;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  if( set_value == NULL ){
    set_value = "__NULL_STR";
  }
  result = add_or_update_hash(parameters->hash, name, set_value);
  result = add_or_update_hash(usages->hash, name, usage);
  result = add_or_update_hash(types->hash, name, "STRING_ARG");
    
  return result;
}

BOOLEAN_T temp_set_mass_type_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 MASS_TYPE_T set_value,  ///< the value to be set -in
 char* usage
 )
{
  BOOLEAN_T result;
  char* value_str;

  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify the value */
  mass_type_to_string(set_value, &value_str);
  
  result = add_or_update_hash(parameters->hash, name, value_str);
  result = add_or_update_hash(usages->hash, name, usage);
  result = add_or_update_hash(types->hash, name, "MASS_TYPE_T");
    
  return result;

}

BOOLEAN_T temp_set_peptide_type_parameter(
 char*     name,  ///< the name of the parameter looking for -in
 PEPTIDE_TYPE_T set_value,  ///< the value to be set -in
 char* usage
  )
{
  BOOLEAN_T result;
  char* value_str;
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  /* stringify the value */
  mass_type_to_string(set_value, &value_str);

  result = add_or_update_hash(parameters->hash, name, value_str);
  result = add_or_update_hash(usages->hash, name, usage);
  result = add_or_update_hash(types->hash, name, "PEPTIDE_TYPE_T");
    
  return result;

}


/**************************************************
 *   OLD SETTERS (public)
 **************************************************
 */


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
    
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }

  // only check if parameter file has already been parsed
  if(parameter_parsed || parameter_initialized){
    if(set_value){
      return update_hash_value(parameters->hash, name, "TRUE");
    }
    else{
      return update_hash_value(parameters->hash, name, "FALSE");
    }
  }

  // if it doesn't already exist(wasn't in the parameter file), add to parameter list
  if(set_value){
    result = add_parameter(name, "TRUE");
  }
  else{
    result = add_parameter(name, "FALSE");
  }
  
  return result;
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
  
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // only check if parameter file has already been parsed
  if(parameter_parsed  || parameter_initialized){
    snprintf(buffer, PARAMETER_LENGTH, "%d", set_value);
    return update_hash_value(parameters->hash, name, buffer);
  }
  
  // if it doesn't already exist(wasn't in the parameter file), add to parameter list
  snprintf(buffer, PARAMETER_LENGTH, "%d", set_value);  
  result = add_parameter(name, buffer);
  
  return result;
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
  
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // convert to string
  sprintf(buffer, "%f", set_value);

  // only check if parameter file has already been parsed
  if(parameter_parsed || parameter_initialized){
    return update_hash_value(parameters->hash, name, buffer);    
  }

  // if it doesn't already exist(wasn't in the parameter file), add to parameter list
  result = add_parameter(name, buffer);
    
  return result;
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
  
  // check if parameters cah be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // only check if parameter file has already been parsed
  if(parameter_parsed  || parameter_initialized){
    return update_hash_value(parameters->hash, name, set_value);
  }
  
  // if it doesn't already exist(wasn't in the parameter file), add to parameter list
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
  // check if parameters has been initlialized
  if(!parameter_initialized && !parameter_parsed){
    carp(CARP_ERROR, "Must initialize parameters before copying");
    return FALSE;
  }
  
  // check if parameters can be changed
  if(!parameter_plasticity){
    carp(CARP_ERROR, "Can't change parameters once they are confirmed");
    return FALSE;
  }
  
  // for required options, there are not in the parameter list, thus must add
  if(required){
    if(add_parameter(name, set_value)){
      return TRUE;
    }
    return FALSE;
  }

  // if exist ovewrite it!
  return update_hash_value(parameters->hash, name, set_value);  
} 

//move to getters section
/**
 * Routines that return crux enumerated types. 
 */
/*
static char* peptide_type_strings[NUMBER_PEPTIDE_TYPES] = {
  "tryptic", "partial", "n-tryptic", "c-tryptic", "not-tryptic", "all"
}; 
*/
BOOLEAN_T string_to_param_type(char* name, PARAMETER_TYPE_T* result ){
  BOOLEAN_T success = TRUE;
  int param_type = convert_enum_type_str(
	 name, -10, parameter_type_strings, NUMBER_PARAMETER_TYPES);
  (*result) = (PARAMETER_TYPE_T)param_type;

  if( param_type == -10 ){
    success = FALSE;
  }
  return success;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

