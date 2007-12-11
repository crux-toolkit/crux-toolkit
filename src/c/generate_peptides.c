/*****************************************************************************
 * \file generate_peptides
 * AUTHOR: Chris Park
 * CREATE DATE: July 17 2006
 * DESCRIPTION: Given a protein fasta sequence database as input, generate a list of peptides in 
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
#include "carp.h"
#include "peptide.h"
#include "peptide_src.h"
#include "protein.h"
#include "database.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "index.h"
#include "generate_peptides_iterator.h"

#define NUM_GEN_PEP_OPTIONS 14
#define NUM_GEN_PEP_ARGS 1

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg){
  char* usage = parse_arguments_get_usage("generate_peptides");
  carp(CARP_FATAL, "incorrect argument %s", arg);
  fprintf(stderr, "%s", usage);
  free(usage);
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
/* Private function declarations */
void print_header();

int main(int argc, char** argv){

  /* Declarations */
  int verbosity;
  /*  int min_length;
  int max_length;
  double min_mass;
  double max_mass;
  MASS_TYPE_T isotopic_mass;
  PEPTIDE_TYPE_T cleavages;
  int missed_cleavages;
  */

  //TODO change ints to bools
  //int use_index;
  BOOLEAN_T output_sequence;
  int trypticity_opt;  //change to print_trypticity
  BOOLEAN_T print_trypticity;
  //char* sort = NULL;   //change to SORT_TYPE_T
  //char* in_file = NULL;
  //create outfile name
  
  long total_peptides = 0;
  GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator = NULL; 
  PEPTIDE_T* peptide = NULL;
    

  /* Set default values for any options here */

  //  int flag_opt = FALSE;
  //  int trypticity_opt = FALSE;
  //  double min_mass = 200;
  //  double max_mass = 7200;
  //  int min_length = 6;
  //  int max_length = 50;
  //  char* cleavages = "tryptic"; 
  //  char* isotopic_mass = "average" ;
  //  int  verbosity = CARP_INFO;
  //  char* use_index = "F";
  //  char* parameter_file = "crux.params";  //remove

  //  int missed_cleavages = FALSE;
  //  char* sort = "none";      // mass, length, lexical, none  
  //  char * in_file = NULL;
  //  char * error_message;
  //  int result = 0;
    
  /* Define optional command line arguments */ 
  int num_options = NUM_GEN_PEP_OPTIONS;
  char* option_list[NUM_GEN_PEP_OPTIONS] = {
    "verbosity",
    "parameter-file",
    "min-length",
    "max-length",
    "min-mass",
    "max-mass",
    "isotopic-mass",
    "cleavages",
    "missed-cleavages",
    "unique-peptides",
    "use-index",
    "output-sequence",
    "output-trypticity",
    "sort"
  };

  /* Define required command-line arguments */
  int num_arguments = NUM_GEN_PEP_ARGS;
  char* argument_list[NUM_GEN_PEP_ARGS] = { "protein input" };

  //TODO make this a debug flag
  set_verbosity_level(CARP_DETAILED_DEBUG);

  /* Prepare parameter.c to read command line, set default option values */
  initialize_parameters();

  /* Set optional and required command-line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments );

  /*
  parse_arguments_set_opt(
    "min-mass", 
    "The minimum neutral mass of the peptides to output.", 
    (void *) &min_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "max-mass", 
    "The maximum neutral mass of the peptides to output.", 
    (void *) &max_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "min-length", 
    "The minimum length of the peptides to output.",
    (void *) &min_length, 
    INT_ARG);

  parse_arguments_set_opt(
    "max-length", 
    "The maximum length of the peptides to output. maximum limit = 255.",
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
    "sort", 
    "Specify the order in which peptides are printed to standard output. none|mass|length|lexical.", 
    (void *) &sort, 
    STRING_ARG);

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
    "use-index", 
    "Specify whether a pre-computed on-disk index should be used for retrieving the peptides. T|F",
    (void *) &use_index, 
    STRING_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);

  parse_arguments_set_opt(
    "output-sequence", 
    "Output the peptide sequence as well as the protein id and start and stop.", 
    (void *) &flag_opt, 
    FLAG_ARG);

  parse_arguments_set_opt(
    "output-trypticity", 
    "Output the peptide trypticiy in output.", 
    (void *) &trypticity_opt, 
    FLAG_ARG);
  */
  /* Define required command line arguments */
  /*
  parse_arguments_set_req(
    "fasta-file", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);
  */

  /* Parse the command line, including optional params file
     includes syntax, type, and bounds checks and dies on error */
  parse_cmd_line_into_params_hash(argc, argv);

  /* Set verbosity */
  verbosity = get_int_parameter("verbosity");
  //why is this segfaulting????
  set_verbosity_level(verbosity);

  /* Get parameter values */
  print_trypticity = get_boolean_parameter("output-trypticity");
  //replace with above
  trypticity_opt = (int)print_trypticity;
  output_sequence = get_boolean_parameter("output-sequence");

  /* Parse the command line */
  /*  if (parse_arguments(argc, argv, 0)) {
    long total_peptides = 0;
    GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator = NULL; 
    PEPTIDE_T* peptide = NULL;
    BOOLEAN_T output_sequence;
    
    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity");
    }

    // first, parse paramter file
    // Next, updates the parameter files with command line options
    parse_update_parameters(parameter_file);
        
    // parameters are now confirmed, can't be changed
    parameters_confirmed();
  */

  print_header();

    // print header line
  /*    printf("# PROTEIN DATABASE: %s\n", get_string_parameter_pointer("fasta-file"));
    printf("# OPTIONS:\n");
    printf("#\tmin-mass: %.2f\n", get_double_parameter("min-mass"));
    printf("#\tmax-mass: %.2f\n", get_double_parameter("max-mass"));
    printf("#\tmin-length: %d\n", get_int_parameter("min-length"));
    printf("#\tmax-length: %d\n", get_int_parameter("max-length"));
    printf("#\tcleavages: %s\n", get_string_parameter_pointer("cleavages"));
    printf("#\tallow missed-cleavages: %s\n", get_string_parameter_pointer("missed-cleavages"));
    printf("#\tsort: %s\n",  get_string_parameter_pointer("sort"));
    printf("#\tisotopic mass type: %s\n", get_string_parameter_pointer("isotopic-mass"));
    printf("#\tverbosity: %d\n", verbosity);
    printf("#\tuse index: %s\n", get_string_parameter_pointer("use-index"));
  */
    
    // create peptide interator
    peptide_iterator = new_generate_peptides_iterator();

  carp(CARP_DETAILED_DEBUG, "Got to here");

    
    // should I output sequence? //ABOVE
    //    output_sequence = get_boolean_parameter("output-sequence");

    // iterate over all peptides
    while(generate_peptides_iterator_has_next(peptide_iterator)){
      ++total_peptides;
      peptide = generate_peptides_iterator_next(peptide_iterator);
      print_peptide_in_format(peptide, output_sequence, 
			      trypticity_opt, stdout);
      
      // free peptide
      free_peptide(peptide);

      if(total_peptides% 10000 == 0){
        carp(CARP_INFO, "reached peptide: %d", total_peptides);
      }
    }
    free_generate_peptides_iterator(peptide_iterator);
    
    // debug purpose
    carp(CARP_INFO, "total peptides: %d", total_peptides);
    free_parameters();
    /*
    exit(0);
  } 
  else {
    char* usage = parse_arguments_get_usage("generate_peptides");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
    exit(1);
  }
    */

  /* successfull exit message */
  carp(CARP_INFO, "crux-create-index finished.");

  exit(0);
}

void print_header(){
    printf("# PROTEIN DATABASE: %s\n", 
	               get_string_parameter_pointer("protein input"));
    printf("# OPTIONS:\n");
    printf("#\tmin-mass: %.2f\n", get_double_parameter("min-mass"));
    printf("#\tmax-mass: %.2f\n", get_double_parameter("max-mass"));
    printf("#\tmin-length: %d\n", get_int_parameter("min-length"));
    printf("#\tmax-length: %d\n", get_int_parameter("max-length"));
    printf("#\tcleavages: %s\n", get_string_parameter_pointer("cleavages"));
    printf("#\tallow missed-cleavages: %s\n", 
                          get_string_parameter_pointer("missed-cleavages"));
    printf("#\tsort: %s\n",  get_string_parameter_pointer("sort"));
    printf("#\tisotopic mass type: %s\n", 
                             get_string_parameter_pointer("isotopic-mass"));
    printf("#\tverbosity: %d\n", get_verbosity_level());
    printf("#\tuse index: %s\n", get_string_parameter_pointer("use-index"));
    
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
