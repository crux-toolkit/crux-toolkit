/*****************************************************************************
 * \file create_protein_index
 * AUTHOR: Chris Park
 * CREATE DATE: September 20 2006
 * DESCRIPTION: Given a protein fasta sequence database as input, generate a protein index file
 *              that contain list of proteins with their protein index and file offset from the
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

int main(int argc, char** argv){
  char * in_file = NULL;
  const char * error_message;
  int result = 0;
  int  verbosity = CARP_INFO;

  /* Define required command line arguments */
  parse_arguments_set_req(
    "protein input filename", 
    "The name of the file (in fasta format) from which to parse proteins.", 
    (void *) &in_file, STRING_ARG);
  
  
  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    set_verbosity_level(verbosity);

    //create protein index if not already present
    if(protein_index_on_disk(in_file)){
      carp(CARP_INFO, "protein index already exist on disk");
    }
    else{
      if(!create_protein_index(in_file)){
        carp(CARP_FATAL, "failed to create protein index on disk");
        exit(1);
      }
    }      
    exit(0);
  } 
  else {
    char* usage = parse_arguments_get_usage("create_protein_index");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
  }
  exit(0);
}
