/*************************************************************************//**
 * \file create_protein_index
 * AUTHOR: Chris Park
 * CREATE DATE: September 20 2006
 * DESCRIPTION: Given a protein fasta sequence database as input,
 * generate a protein index file 
 *              that contain list of proteins with their protein index
 *              and file offset from the 
 *              the fasta file
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
#include "parse_arguments.h"
#include "index.h"
#include "protein_index.h"

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg){
  char* usage = parse_arguments_get_usage("create_peptide_index");
  carp(CARP_FATAL, "incorrect argument %s\n%s", arg, usage);
}

int main(int argc, char** argv){
  char * in_file = NULL;
  char* type_of_file = "binary";
  char * error_message;
  int result = 0;
  int  verbosity = CARP_INFO;
  BOOLEAN_T is_binary = TRUE;

  /* Define optional commands */
  parse_arguments_set_opt(
    "output-type", 
    "the type of protein file to create, binary|index", 
    (void *) &type_of_file, STRING_ARG);


  /* Define required command line arguments */
  parse_arguments_set_req(
    "protein input filename", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);
  
  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    set_verbosity_level(verbosity);
    
    // check which file to create
    if(strcmp(type_of_file, "binary") == 0){
      is_binary = TRUE;
    }
    else if(strcmp(type_of_file, "index") == 0){
      is_binary = FALSE;
    }
    else{
      wrong_command(type_of_file);
    }
    
    // create protein index if not already present
    if(protein_index_on_disk(in_file, is_binary)){
      carp(CARP_INFO, "protein %s fasta file already exist on disk", type_of_file);
      exit(0);
    }
    
    // shoudl I create the binary fasta file?
    if(is_binary){
      if(!create_binary_fasta(in_file)){
        carp(CARP_FATAL, "failed to create binary fasta file on disk");
      }
    }
    else{// create protein index file
      if(!create_protein_index(in_file)){
        carp(CARP_FATAL, "failed to create protein index on disk");
      }
    }    
    exit(0);
  } 
  else {
    char* usage = parse_arguments_get_usage("create_protein_index");
    result = parse_arguments_get_error(&error_message);
    carp(
      CARP_FATAL, 
      "Error in command line. Error # %d\nMessage %s\n%s", 
      result,
      error_message,
      usage
    );
  }
  exit(0);
}
