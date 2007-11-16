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
//#include "parse_arguments.h"
#include "index.h"
#include "protein_index.h"
#include "parameter.h"

#define NUM_INDEX_OPTIONS 9
#define NUM_INDEX_ARGS 1

/**
 * when wrong command is seen carp, and exit
 */
//TODO change this to bad_param_value
//     and use for testing string_to_type params
/*
void wrong_command(char* arg){
  //get_usage should now already "know" what is legal
  char* usage = parse_arguments_get_usage("create_index");
  carp(CARP_FATAL, "incorrect argument %s", arg);
  fprintf(stderr, "%s", usage);
  free(usage); //??
  exit(1);
}
*/
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
  int verbosity;
  double min_mass;
  double max_mass;
  int min_length;
  int max_length;
  int missed_cleavages; 
  PEPTIDE_TYPE_T peptide_type;
  //where is the unique vs redundant variable?
  MASS_TYPE_T mass_type;

  double mass_range;
  PEPTIDE_CONSTRAINT_T* constraint;
  char* in_file = NULL;
  INDEX_T* crux_index;
  char* binary_fasta_file = NULL;

  //to be deleted
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
    //    "unique-peptides",
    "missed-cleavages"//,
    //    "mass-range"
  };
  /*BF Define required command line arguments */ 
  //TODO add to this index name
  int num_arguments = NUM_INDEX_ARGS;
  char* argument_list[NUM_INDEX_ARGS] = { "protein fasta file" };


  /* For debugging of parameter processing */
  /* TODO make this dependant on a compile flag */
  set_verbosity_level(CARP_DETAILED_DEBUG);  
  carp(CARP_DETAILED_DEBUG, "Starting create_index");

  // connect various signals to our clean-up function
  signal( SIGTERM, clean_up );
  signal( SIGINT, clean_up );
  signal( SIGQUIT, clean_up );
  signal( SIGHUP, clean_up ); 

  
  /* set up parameters and their defaults in parameter.c */
  initialize_parameters();

  /* Define optional and required command line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  /* Parse the command line, including the optional params file */
  /* does sytnax, type, bounds checking and dies if neccessessary */
  parse_cmd_line_into_params_hash(argc, argv);

  //move this to parameter.c
  // set verbosity
  verbosity = get_int_parameter("verbosity");
  if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
    set_verbosity_level(verbosity);
  }
  else{
    // wrong_command("verbosity");
  }
  
  //  carp(CARP_DETAILED_DEBUG, "Did we get to here?");
  
  /* Get parameter values */
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
  //mass_range = get_double_parameter("mass-range"); //int?
  //check that range > 0
  //change input to take int 0-max-length
  missed_cleavages = get_boolean_parameter("missed-cleavages");
  
  
  // FIXME may add additional types such as non-tryptic or partially-tryptic
  peptide_type = get_peptide_type_parameter("cleavages");
  
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
    //wrong_command(mass_type_str);
  }
  
  // peptide constraint
  constraint = new_peptide_constraint(peptide_type, min_mass, max_mass, 
				      min_length, max_length, 
				      missed_cleavages, mass_type);
  
  // check if input file exist
  in_file = get_string_parameter("protein fasta file");
  carp(CARP_DETAILED_DEBUG,"Input file name is '%s'\n", in_file);
  if(access(in_file, F_OK)){
    carp(CARP_FATAL, "The file \"%s\" does not exist " 
	 "(or is not readable or is empty).", in_file);
    exit(1);
  }
  
  
  mass_range = (max_mass - min_mass)/MAX_INDEX_FILES;
  
  //FINALLY! the computation begins
  
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
}
