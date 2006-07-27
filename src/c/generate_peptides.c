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
#include "protein_peptide_association.h"
#include "protein.h"
#include "database.h"
#include "parse_arguments.h"

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

int main(int argc, char** argv){

  /* Set default values for any options here */
  int flag_opt = FALSE;
  double min_mass = 200;
  double max_mass = 2400;
  int min_length = 6;
  int max_length = 50;
  char* cleavages = "tryptic"; 
  char* isotopic_mass = "average" ;
  int  verbosity = CARP_MAX;
  char* redundancy = "redundant";

  MASS_TYPE_T mass_type = AVERAGE;
  PEPTIDE_TYPE_T peptide_type = TRYPTIC;
  int missed_cleavages = FALSE;
  char* sort = "none";      // mass or length or none
  char * in_file = NULL;
  const char * error_message;
  int result = 0;
  BOOLEAN_T is_unique = FALSE;

  /* Define optional command line arguments */ 
  parse_arguments_set_opt(
    "output-sequence", 
    "Output the peptide sequence as well as the protein id and start and stop.", 
    (void *) &flag_opt, 
    FLAG_ARG);

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
    "The maximum length of the peptides to output.",
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
    "Specify the order in which peptides are printed to standard output. ", 
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
    "redundancy", 
    "Specify whether peptides that come from different proteins yet with identical sequences should appear on separate lines or on the same line. redundant|unique.",
    (void *) &redundancy, 
    STRING_ARG);

  /* Define required command line arguments */
  parse_arguments_set_req(
    "protein input filename", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);


  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    PEPTIDE_CONSTRAINT_T* constraint;
    DATABASE_PEPTIDE_ITERATOR_T* iterator;
    DATABASE_T* database;
    PEPTIDE_T* peptide;
   
    //FIXME may add additional types such as non-trypticc or partially-tryptic
    if(strcmp(cleavages, "all")==0){
      peptide_type = ANY_TRYPTIC;
    }
    else if(strcmp(cleavages, "tryptic")==0){
      peptide_type = TRYPTIC;
    }
    else if(strcmp(cleavages, "partial")==0){
      peptide_type = PARTIALLY_TRYPTIC;
    }
    else{
      wrong_command(cleavages);
    }
    
    //determine isotopic mass option
    if(strcmp(isotopic_mass, "average")==0){
      mass_type = AVERAGE;
    }
    else if(strcmp(isotopic_mass, "mono")==0){
      mass_type = MONO;
    }
    else{
      wrong_command(isotopic_mass);
    }
   
    //determine redundancy option
    if(strcmp(redundancy, "redundant")==0){
      is_unique = FALSE;
    }
    else if(strcmp(redundancy, "unique")==0){
      is_unique = TRUE;
    }
    else{
      wrong_command(redundancy);
    }

    //set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity");
    }
  
    //peptide constraint
    constraint = new_peptide_constraint(peptide_type, min_mass, max_mass, min_length, max_length, missed_cleavages, mass_type);
 
    //check if input file exist
    if(access(in_file, F_OK)){
      carp(CARP_FATAL, "The file \"%s\" does not exist (or is not readable, or is empty).", in_file);
      exit(1);
    }
    
    //print header line
    printf("# PROTEIN DATABASE: %s\n", in_file);
    printf("# OPTIONS:\n");
    printf("#\tmin-mass: %.2f\n", min_mass);
    printf("#\tmax-mass: %.2f\n", max_mass);
    printf("#\tmin-length: %d\n", min_length);
    printf("#\tmax-length: %d\n", max_length);
    printf("#\tcleavages: %s\n", cleavages);
    if(missed_cleavages){
      printf("#\tallow missed-cleavages: TRUE\n");
    }
    else{
      printf("#\tallow missed-cleavages: FALSE\n");
    }
    printf("#\tsort: %s (needs to be implemented)\n", sort);
    if(mass_type == AVERAGE){
      printf("#\tisotopic mass type: average\n");
    }
    else{
      printf("#\tisotopic mass type: mono\n");
    }
    printf("#\tverbosity: %d\n", verbosity);
    if(is_unique){
      printf("#\tredundancy: unique\n");
    }
    else{
      printf("#\tredundancy: redundant\n");
    }

    //create a new database
    database = new_database(in_file);
    //create peptide iterator
    iterator = new_database_peptide_iterator(database, constraint);
    
    //check if any peptides are found
    if(!database_peptide_iterator_has_next(iterator)){
      carp(CARP_WARNING, "no matches found\n");
    }
    else{
      int progress = 0;
      fprintf(stderr, "Running:\n");
      //print each peptide
      while(database_peptide_iterator_has_next(iterator)){
        show_progress(&progress);
        peptide = database_peptide_iterator_next(iterator);
        print_peptide_in_format(peptide, flag_opt, stdout);
        free_peptide(peptide);
      }
      putc('\n', stderr);
      fflush(stderr);
    }
      //free database, iterator, constraint
      free_database_peptide_iterator(iterator);
      free_peptide_constraint(constraint);
      free_database(database);      
      exit(0);
  } 
  else {
    char* usage = parse_arguments_get_usage("generate_peptides");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
  }
  exit(0);
}
