/*****************************************************************************
 * \file create_index
 * AUTHOR: Chris Park
 * CREATE DATE: August 10 2006
 * DESCRIPTION: Given a protein fasta sequence database as input, generate crux_index files
 *              that contain list of peptides in 
 *              the database that meet certain criteria (e.g. mass, length, trypticity) as output.
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include <signal.h>
#include "carp.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
//all of parse_arguments should be done through parameter.c
#include "parse_arguments.h"
#include "index.h"
#include "protein_index.h"
#include "parameter.h"

#define NUM_INDEX_OPTIONS 11
#define NUM_INDEX_ARGS 1

/**
 * when wrong command is seen carp, and exit
 */
//TODO change this to bad_param_value
//     and use for testing string_to_type params
void wrong_command(char* arg){
  //get_usage should now already "know" what is legal
  char* usage = parse_arguments_get_usage("create_index");
  carp(CARP_FATAL, "incorrect argument %s", arg);
  fprintf(stderr, "%s", usage);
  free(usage); //??
  exit(1);
}

/**
 * progress indicator, displays a spinning | to standout
 */
void show_progress(int* num){
  putc('\b', stderr);
  if(*num / 150 == 37){
    putc('|', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 74){
    putc('/', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 111){
    putc('-', stderr);
    fflush(stderr);
  }
  else if(*num / 150 == 149){
    putc('\\', stderr);
    fflush(stderr);
    *num = 0;
    return;
  }
  ++*num;
}
      
int main(int argc, char** argv){

  /* Declarations */
  //TODO put these in the decided order
  int verbosity, min_length, max_length, missed_cleavages; 
  double min_mass, max_mass;
  PEPTIDE_TYPE_T peptide_type;
  MASS_TYPE_T mass_type;
  PEPTIDE_CONSTRAINT_T* constraint;
  INDEX_T* crux_index;
  char * in_file = NULL;
  char* binary_fasta_file = NULL;

  //to be deleted
  double mass_range;
  char* mass_type_str;

  /* Define optional command line arguments */ 

  int num_options = NUM_INDEX_OPTIONS;
  char* option_list[NUM_INDEX_OPTIONS] = { 
    "verbosity",
    "parameter-file", 
    "min-mass", 
    "max-mass", 
    "min-length", 
    "max-length", 
    "cleavages", 
    "isotopic-mass",
    "unique-peptides",
    "missed-cleavages",
    "mass-range"
  };
  /*BF Define required command line arguments */ 
  //TODO add to this index name
  int num_arguments = NUM_INDEX_ARGS;
  char* argument_list[NUM_INDEX_ARGS] = { "protein fasta file" };


  /* For debugging of parameter processing */
  /* make this dependant on a compile flag */
  set_verbosity_level(CARP_DETAILED_DEBUG);  
  carp(CARP_DETAILED_DEBUG, "Starting create_index");

  // connect various signals to our clean-up function
  signal( SIGTERM, clean_up );
  signal( SIGINT, clean_up );
  signal( SIGQUIT, clean_up );
  signal( SIGHUP, clean_up ); 

  
  /* Set default values for any options here */
  //WAIT: set these after parsing command line
  //set these all in the form
  //double min_mass = get_double_parameter("min-mass");
  /*  double min_mass = 200;
  double max_mass = 7200;
  double mass_range = 1;  
  int min_length = 6;
  int max_length = 50;
  char* cleavages = "tryptic"; 
  char* isotopic_mass = "average" ;
  int  verbosity = CARP_INFO;
  char* binary_fasta_file = NULL;
  */  //  char* parameter_file = "crux.params";
  /*  char* parameter_file = NULL;
  
  MASS_TYPE_T mass_type = AVERAGE;
  PEPTIDE_TYPE_T peptide_type = TRYPTIC;
  int missed_cleavages = FALSE;
  char * in_file = NULL;
*/
  /*const*/ //char * error_message;
  //int result = 0;
  
  /* set up parameters and their defaults in parameter.c */
  initialize_parameters();

  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  //at what level should the error checking be done?
  //deepest: parse_arguments -- dies with message about the bad option
  //mid:     parameter -- can catch parse_argument messages, could check 
  //                      extra types and bounds (if same for all exe's),
  //                      could use parse_argument to generate usage
  //here:    exe -- would require that parameter catch and keep all parse_arg
  //                error info, each exe has to check same bounds
  parse_cmd_line_into_params_hash(argc, argv);

  //NEXT: get values back...move all bounds checking to params
  /*
  parse_arguments_set_opt(
    "mass-range", 
    "The mass range contained in each index file.",
    (void *) &mass_range, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "min-mass", 
    "The minimum mass of the peptides to put in index file.", 
    (void *) &min_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "max-mass", 
    "The maximum mass of the peptides to output in index file.", 
    (void *) &max_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "min-length", 
    "The minimum length of the peptides to output in index file.",
    (void *) &min_length, 
    INT_ARG);

  parse_arguments_set_opt(
    "max-length", 
    "The maximum length of the peptides to output in index file. maximum limit = 255.",
    (void *) &max_length, 
    INT_ARG);

  parse_arguments_set_opt(
    "cleavages", 
    "Type of cleavages to allow. tryptic|partial|all.", 
    (void *) &cleavages, 
    STRING_ARG);

  parse_arguments_set_opt(
    "missed-cleavages", 
    "Allow missed cleavage sites with in a peptide. ",
    (void *) &missed_cleavages, 
    FLAG_ARG);
 
  parse_arguments_set_opt(
    "isotopic-mass", 
    "Specify the type of isotopic masses to use when calculating the peptide mass. average|mono.",
    (void *) &isotopic_mass, 
    STRING_ARG);

  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG); 
  */
  /* Define required command line arguments */

  //now as 
  //int num_args = 1;
  //char** arg_list = ("protein input filename");
  //select_cmd_line_arguments( arg_list, num_args );
  /*  parse_arguments_set_req(
    "protein input filename", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);
  */

  /* Parse the command line */

  //now params.c parses command line
  //if( ! parse_arguments_into_params( argc, argv) ){
  //   char* usage = get_usage();
  //   ... print error and die
  //}

  //  if (parse_arguments(argc, argv, 0)) {
  //move these?
  //PEPTIDE_CONSTRAINT_T* constraint;
  //INDEX_T* crux_index;

  //char* binary_fasta_file = NULL;
    //char * in_file = NULL;
  /*const*/// char * error_message;
  //int result = 0;

  
  //move this to parameter.c
    // set verbosity
    verbosity = get_int_parameter("verbosity");
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity");
    }
    
    //  carp(CARP_DETAILED_DEBUG, "Did we get to here?");
  
    // get parameters
    // move all type-checking and bounds checking to parameter.c
    // note: parse_argument only knows native types
    //       parameter knows all types ahd has get_ calls for each
    //       parameter uses string_to_type functions (but they are also available to others)
    min_mass = get_double_parameter("min-mass");
    max_mass = get_double_parameter("max-mass");
    //check that min > 0 and min < max and max < 255
    min_length = get_int_parameter("min-length");
    max_length = get_int_parameter("max-length");
    //check that min > 0 and min < max and max < 255
    mass_range = get_double_parameter("mass-range"); //int?
    //check that range > 0
    //change input to take int 0-max-length
    missed_cleavages = get_boolean_parameter("missed-cleavages");

    //already done
    // parse and update parameters
    /*    parse_update_parameters(parameter_file);

    // parameters are now confirmed, can't be changed
    parameters_confirmed();
    */
    /******* All parameters must be taken through get_*_parameter() method ******/

    // FIXME may add additional types such as non-tryptic or partially-tryptic
    peptide_type = get_peptide_type_parameter("cleavages");
    //or rather than parameter, make a string_to_peptide_type(char*, PEPTIDE_TYPE*) function
    /*
    if(strcmp(get_string_parameter_pointer("cleavages"), "all")==0){
      peptide_type = ANY_TRYPTIC;
    }
    else if(strcmp(get_string_parameter_pointer("cleavages"), "tryptic")==0){
      peptide_type = TRYPTIC;
    }
    else if(strcmp(get_string_parameter_pointer("cleavages"), "partial")==0){
      peptide_type = PARTIALLY_TRYPTIC;
      }
    else{
      wrong_command(cleavages);
    }
    */

    //again, move bounds checking
    // check if maximum length is with in range <= 255
    //put above
    if(max_length > 255){
      carp(CARP_FATAL, "Maximum length:%d over limit of 255.", max_length);
      exit(1);
    }
    
    // mass_range = get_double_parameter("mass-range");
    /*if(compare_float(mass_range, 0) == 0){
      carp(CARP_FATAL, "mass_range:%d must be greater than 0.", mass_range);
      exit(1);
    }*/
    
    // determine isotopic mass option
    // write a BOOLEAN_T string_to_mass_type(char*, MASS_TYPE_T*) function
    mass_type_str = get_string_parameter("isotopic-mass");
    //MASS_TYPE_T mass_type;
    if(strcmp(mass_type_str, "average")==0){
      mass_type = AVERAGE;
    }
    else if(strcmp(mass_type_str, "mono")==0){
      mass_type = MONO;
    }
    else{
      //      wrong_command(isotopic_mass);
      wrong_command(mass_type_str);
    }

    // peptide constraint
    constraint = new_peptide_constraint(peptide_type, min_mass, max_mass, 
        min_length, max_length, missed_cleavages, mass_type);
    
    // check if input file exist
    in_file = get_string_parameter("protein fasta file");
    carp(CARP_DETAILED_DEBUG,"Input file name is '%s'\n", in_file);
    if(access(in_file, F_OK)){
      carp(CARP_FATAL, "The file \"%s\" does not exist " 
          "(or is not readable or is empty).", in_file);
      exit(1);
    }
    
    //FINALLY! the computation begins
    mass_range = (max_mass - min_mass)/MAX_INDEX_FILES;

    // create new index object
    crux_index = 
      new_index(in_file,
                constraint,
                mass_range
                );
    
    // create crux_index files
    if(!create_index(crux_index)){
      die("Failed to create index");
    }
    
    // free index(frees constraint together);
    free_index(crux_index);     
    free(binary_fasta_file);
    free_parameters();
    //successfull exit message
    exit(0);
    //  } 
    /*  else {
    char* usage = parse_arguments_get_usage("create_index");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
    exit(1);
    }*/
  exit(0);
}
