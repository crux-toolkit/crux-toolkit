/**
 * \file crux-main.cpp
 */
/*
 AUTHOR: Barbara Frewen
 CREATE DATE: November 24, 2008
 DESCRIPTION: The starting point for what were previously three
         separate protgrams--create-index, search-for-matches,
         analyze-matches.  Usage is
            crux <operation> <options> <arguments>
         where operation is create-index, search, compute-q-values, or
         q-ranker.     
 REVISION: $Revision: 1.2 $
*/

#include "crux-main.h"

const char* usage_str = "Usage: crux <command> [options] <argument>\n"
"Commands:\n"
"  create-index        Create an index for all peptides in a fasta file.\n"
"  search-for-matches  Search a collection of spectra against a sequence\n"
"                      database, returning a collection of peptide-spectrum\n"
"                      matches (PSMs) scored by XCorr.\n"
"  sequest-search      Similar to search-for-matches but use Sp as a \n"
"                      preliminary score followed by XCorr.\n"
"  compute-q-values    Assign a q-value, which is a statistical confidence\n"
"                      measure that accounts for multiple testing, to each\n"
"                      PSM in a given set.\n" 
"  percolator          Analyze a collection of PSMs to target and decoy\n"
"                      sequences using the percolator algorithm.\n"
"  q-ranker            Analyze a collection of PSMs using the Q-ranker\n"
"                      algorithm.\n"
"  print-processed-spectra\n"
"                      Write a new ms2 file with all of the same spectra\n"
"                      with only the peaks used for computing xcorr.\n"
"Options and arguments:\n"
"  Specific to each command. Type 'crux <command>' to get details.\n"
;
int main(int argc, char** argv){

  // check the syntax for crux <operation>
  if( argc < 2 ){
    carp(CARP_FATAL, usage_str);
  }

  // determine the operation
  char* op_string = argv[1];
  COMMAND_T command = string_to_command_type(op_string);

  // call the appropriate function 
  // passing the command line minus the first token ('crux')
  switch(command){
  case INDEX_COMMAND:
    create_index_main(argc-1, argv+1);
    break;

  case SEARCH_COMMAND:
    search_main(argc-1, argv+1);
    break;

  case SEQUEST_COMMAND:
    sequest_search_main(argc-1, argv+1);
    break;

  case QVALUE_COMMAND:
    qvalue_main(argc-1, argv+1);
    break;

  case QRANKER_COMMAND:
    qranker_main(argc-1, argv+1);
    break;

  case PERCOLATOR_COMMAND:
    percolator_main(argc-1, argv+1);
    break;

  case PROCESS_SPEC_COMMAND:
    print_processed_spectra_main(argc-1, argv+1);
    break;
  case INVALID_COMMAND:
    carp(CARP_FATAL, "Invalid command '%s'\n%s", op_string, usage_str);
    break;

  default:
    carp(CARP_FATAL, "Unknown command type.");
    break;

  }

  exit (0);
}// end main















