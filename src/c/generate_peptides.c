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
#include "index.h"

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
  char* redundancy = "unique";
  char* use_index = "F";

  BOOLEAN_T use_index_boolean = FALSE;
  MASS_TYPE_T mass_type = AVERAGE;
  PEPTIDE_TYPE_T peptide_type = TRYPTIC;
  int missed_cleavages = FALSE;
  char* sort = "none";      // mass, length, lexical, none  
  char * in_file = NULL;
  const char * error_message;
  int result = 0;
  BOOLEAN_T is_unique = FALSE;
  SORT_TYPE_T sort_type = NONE;

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
    "redundancy", 
    "Specify whether peptides that come from different proteins yet with identical sequences should appear on separate lines or on the same line. redundant|unique.",
    (void *) &redundancy, 
    STRING_ARG);
  
  parse_arguments_set_opt(
    "use-index", 
    "Specify whether a pre-computed on-disk index should be used for retrieving the peptides. T|F",
    (void *) &use_index, 
    STRING_ARG);


  /* Define required command line arguments */
  parse_arguments_set_req(
    "protein input filename", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);
  

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    PEPTIDE_CONSTRAINT_T* constraint;
    DATABASE_PEPTIDE_ITERATOR_T* iterator = NULL;
    DATABASE_SORTED_PEPTIDE_ITERATOR_T* sorted_iterator = NULL;
    DATABASE_T* database = NULL;
    PEPTIDE_T* peptide = NULL;
    INDEX_T* index = NULL;
    INDEX_PEPTIDE_ITERATOR_T* index_peptide_iterator = NULL;

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
    
    //check if maximum length is with in range <= 255
    if(max_length > 255){
      carp(CARP_FATAL, "maximum length:%d over limit 255.", max_length);
      exit(1);
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

    //determine sort type option
    if(strcmp(sort, "mass")==0){
      sort_type = MASS;
    }
    else if(strcmp(sort, "length")==0){
      sort_type = LENGTH;
    }
    else if(strcmp(sort, "lexical")==0){
      sort_type = LEXICAL;
    }
    else if(strcmp(sort, "none")==0){
      sort_type = NONE;
    }
    else{
      wrong_command(sort);
    }
    
    //determine use index command
    if(strcmp(use_index, "F")==0){
      use_index_boolean = FALSE;
    }
    else if(strcmp(use_index, "T")==0){
      use_index_boolean = TRUE;
    }
    else{
      wrong_command(use_index);
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
    printf("#\tsort: %s\n", sort);
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
    if(use_index_boolean){
      printf("#\tuse index: TRUE\n");
    }
    else{
      printf("#\tuse index: FALSE\n");
    }

    /***********************
     * use index file
     **********************/
    if(use_index_boolean){
      if((sort_type != MASS && sort_type != NONE) || !is_unique){
        carp(CARP_ERROR, " when using index, cannot sort other than by mass and only returns unique peptides");
        carp(CARP_ERROR, "failed to perform search");
        exit(1);
      }

      index = new_search_index(in_file, constraint);

      if(index != NULL){
        //create index peptide interator
        index_peptide_iterator = new_index_peptide_iterator(index);//, ok_seq);
        
        //iterate over all peptides
        while(index_peptide_iterator_has_next(index_peptide_iterator)){
          peptide = index_peptide_iterator_next(index_peptide_iterator);
          print_peptide_in_format(peptide, flag_opt, stdout);
          free_peptide(peptide);
        }
        
        free_index(index);
        free_index_peptide_iterator(index_peptide_iterator);
      }
      else{
        carp(CARP_ERROR, "failed to perform search");
        exit(1);
      }
    }
    /*********************************************
     *read in from fasta file, don't use index file
     ************************************************/
    else{
      //create a new database
      database = new_database(in_file);
      
      //no sort, redundant
      if(!is_unique && sort_type == NONE){ 
        //create peptide iterator
        iterator = new_database_peptide_iterator(database, constraint);
        
        //check if any peptides are found
        if(!database_peptide_iterator_has_next(iterator)){
          carp(CARP_WARNING, "no matches found");
        }
        else{
          //print each peptide
          while(database_peptide_iterator_has_next(iterator)){
            peptide = database_peptide_iterator_next(iterator);
            print_peptide_in_format(peptide, flag_opt, stdout);
            free_peptide(peptide);
          }
        }
        //free iterator
        free_database_peptide_iterator(iterator);
      }      
      //sort or check for unique
      else{
        //only sort, by default will be sorted by mass
        if(sort_type == NONE){
          //create peptide iterator
          sorted_iterator = 
            new_database_sorted_peptide_iterator(database, constraint, MASS, TRUE);
        }
        //create peptide iterator
        else{
          sorted_iterator = 
            new_database_sorted_peptide_iterator(database, constraint, sort_type, is_unique);
        }
        
        //check if any peptides are found
        if(!database_sorted_peptide_iterator_has_next(sorted_iterator)){
          carp(CARP_WARNING, "no matches found");
        }
        else{
          //print each peptide
          while(database_sorted_peptide_iterator_has_next(sorted_iterator)){
            peptide = database_sorted_peptide_iterator_next(sorted_iterator);
            print_peptide_in_format(peptide, flag_opt, stdout);
            free_peptide(peptide);
          }
        }
        //free iterator
        free_database_sorted_peptide_iterator(sorted_iterator);
      }
      //free database, iterator, constraint
      free_peptide_constraint(constraint);
      free_database(database);      
      exit(0);
    }
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
