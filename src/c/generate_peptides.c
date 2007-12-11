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

/* Private function declarations */
void print_header();

int main(int argc, char** argv){

  /* Declarations */
  int verbosity;

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

  print_header();

  // create peptide interator
  peptide_iterator = new_generate_peptides_iterator();
  
  carp(CARP_DETAILED_DEBUG, "Got to here");
  
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
