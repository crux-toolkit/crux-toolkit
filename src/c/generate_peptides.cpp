/*************************************************************************//**
 * \file generate_peptides.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: July 17 2006
 * \brief Given a protein fasta sequence database as input,
 * generate a list of peptides in the database that meet certain
 * criteria (e.g. mass, length, trypticity) as output. 
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
#include "Protein.h"
#include "database.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "index.h"
#include "generate_peptides_iterator.h"

static const int NUM_GEN_PEP_OPTIONS = 15;
static const int NUM_GEN_PEP_ARGS = 1;

/* Private function declarations */
void print_header();

int main(int argc, char** argv){

  /* Declarations */
  //int verbosity;
  BOOLEAN_T output_sequence;
  //  BOOLEAN_T print_trypticity = FALSE;
  BOOLEAN_T use_index;
  char* filename;
  
  long total_peptides = 0;
  //GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator = NULL; 
  MODIFIED_PEPTIDES_ITERATOR_T* peptide_iterator = NULL; 
  PEPTIDE_T* peptide = NULL;
  DATABASE_T* database = NULL;
  INDEX_T* index = NULL;
    
  /* Define optional command line arguments */ 
  int num_options = NUM_GEN_PEP_OPTIONS;
  const char* option_list[NUM_GEN_PEP_OPTIONS] = {
    "version",
    "verbosity",
    "parameter-file",
    "min-length",
    "max-length",
    "min-mass",
    "max-mass",
    "isotopic-mass",
    "enzyme", 
    "custom-enzyme", 
    "digestion", 
    //    "cleavages",
    "missed-cleavages",
    "unique-peptides",
    //"use-index",
    "output-sequence",
    //"output-trypticity",
    "sort"
  };

  /* Define required command-line arguments */
  int num_arguments = NUM_GEN_PEP_ARGS;
  const char* argument_list[NUM_GEN_PEP_ARGS] = { "protein database" };

  //TODO make this a debug flag
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  set_verbosity_level(CARP_ERROR);

  /* Prepare parameter.c to read command line, set default option values */
  initialize_parameters();

  /* Set optional and required command-line arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments );


  /* Parse the command line, including optional params file
     includes syntax, type, and bounds checks and dies on error */
  parse_cmd_line_into_params_hash(argc, argv, "crux-generate-peptides");

  /* Set verbosity */
  //verbosity = get_int_parameter("verbosity");
  //set_verbosity_level(verbosity);

  /* Get parameter values */
  //  print_trypticity = get_boolean_parameter("output-trypticity");
  output_sequence = get_boolean_parameter("output-sequence");
  filename = get_string_parameter("protein database");
  //use_index = get_boolean_parameter("use-index");
  use_index = is_directory(filename);

  if( use_index == TRUE ){
    index = new_index_from_disk(filename);//, 
  }else{
    database = new_database(filename, FALSE); // not memmapped
  }
  free(filename);


  // get list of mods
  PEPTIDE_MOD_T** peptide_mods = NULL;
  int num_peptide_mods = generate_peptide_mod_list( &peptide_mods );
  carp(CARP_DEBUG, "Got %d peptide mods", num_peptide_mods);

  /* Generate peptides and print to stdout */

  print_header(); // TODO: add mods info

  // one iterator for each peptide mod
  int mod_idx = 0;
  for(mod_idx = 0; mod_idx < num_peptide_mods; mod_idx++){
    carp(CARP_DETAILED_DEBUG, "Using peptide mod %d with %d aa mods", 
         mod_idx, peptide_mod_get_num_aa_mods(peptide_mods[mod_idx]));
    // create peptide iterator
    //peptide_iterator = new_generate_peptides_iterator();
    peptide_iterator = new_modified_peptides_iterator(peptide_mods[mod_idx],
                                                      index,
                                                      database);
  
    //print_header(); // TODO: add mods info

    // iterate over all peptides
    int orders_of_magnitude = 1000; // for counting when to print
    //while(generate_peptides_iterator_has_next(peptide_iterator)){
    while(modified_peptides_iterator_has_next(peptide_iterator)){
      ++total_peptides;
      //peptide = generate_peptides_iterator_next(peptide_iterator);
      peptide = modified_peptides_iterator_next(peptide_iterator);
      print_peptide_in_format(peptide, output_sequence, 
                              //print_trypticity, stdout);
                               stdout);
    
      // free peptide
      free_peptide(peptide);
    
      if(total_peptides % orders_of_magnitude == 0){
        if( (total_peptides)/10 == orders_of_magnitude){
          orders_of_magnitude *= 10;
        }
        carp(CARP_INFO, "Reached peptide %d", total_peptides);
      }
    }// last peptide
    //free_generate_peptides_iterator(peptide_iterator);
    free_modified_peptides_iterator(peptide_iterator);

  }// last peptide modification

  // debug purpose
  carp(CARP_INFO, "total peptides: %d", total_peptides);
  free_parameters();

  /* successfull exit message */
  carp(CARP_INFO, "crux-generate-peptides finished.");

  exit(0);
}

void print_header(){
  BOOLEAN_T bool_val;

  char* database_name = get_string_parameter("protein database");
  printf("# PROTEIN DATABASE: %s\n", database_name);

  printf("# OPTIONS:\n");
  printf("#\tmin-mass: %.2f\n", get_double_parameter("min-mass"));
  printf("#\tmax-mass: %.2f\n", get_double_parameter("max-mass"));
  printf("#\tmin-length: %d\n", get_int_parameter("min-length"));
  printf("#\tmax-length: %d\n", get_int_parameter("max-length"));
  printf("#\tenzyme: %s\n", get_string_parameter_pointer("enzyme"));
  printf("#\tdigestion: %s\n", get_string_parameter_pointer("digestion"));
  //printf("#\tcleavages: %s\n", get_string_parameter_pointer("cleavages"));
  
  bool_val = get_boolean_parameter("missed-cleavages");
  printf("#\tallow missed-cleavages: %s\n", boolean_to_string(bool_val));
  printf("#\tsort: %s\n",  get_string_parameter_pointer("sort"));
  printf("#\tisotopic mass type: %s\n", 
         get_string_parameter_pointer("isotopic-mass"));
  printf("#\tverbosity: %d\n", get_verbosity_level());

  //bool_val = get_boolean_parameter("use-index");
  bool_val = is_directory(database_name);
  printf("#\tuse index: %s\n", boolean_to_string(bool_val));
  //get_string_parameter_pointer("use-index"));
  
  AA_MOD_T** aa_mod_list = NULL;
  int num_aa_mods = get_all_aa_mod_list(&aa_mod_list);
  int mod_idx = 0;
  for(mod_idx=0; mod_idx < num_aa_mods; mod_idx++){
    printf("#\tmodification: ");
    char* mod_str = aa_mod_to_string(aa_mod_list[mod_idx]);
    printf("%s\n", mod_str);
    free(mod_str);
  }
  free(database_name);
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
