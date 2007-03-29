/*****************************************************************************
 * \file search_spectrum
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * DESCRIPTION: Given as input an ms2 file, a spectrum scan number within that file, 
 * a sequence database, and an optional parameter file, 
 * search the spectrum against the peptides in the sequence database, and return high scoring peptides. 
 * REVISION: 
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "carp.h"
#include "peptide.h"
#include "protein.h"
#include "parse_arguments.h"
#include "parameter.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "crux-utils.h"
#include "scorer.h"
#include "objects.h"
#include "parameter.h"
#include "match.h"
#include "match_collection.h"
/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){
  char* usage = parse_arguments_get_usage("search_spectrum");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  //print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  //FIXME uncomment this print if want to print usage whenever error message is printed
  //fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

int main(int argc, char** argv){
  /* Set default values for any options here */

  //optional
  int charge = 2;
  char* perlim_score_type = "sp";
  char* score_type = "sp";
  char* parameter_file = NULL;
  int verbosity = CARP_ERROR;
  
  //required
  char* ms2_file = NULL;
  int scan_num = 0;
  char* fasta_file = NULL;
  double mass_window = 3;

  //parsing variables
  int result = 0;
  char * error_message;

  /* Define optional command line arguments */ 
  parse_arguments_set_opt(
    "charge", 
    "The peptide charge. 1|2|3",
    (void *) &charge, 
    INT_ARG);
  
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
  
  parse_arguments_set_opt(
    "score-type", 
    "The type of scoring function to use. sp | xcorr",
    (void *) &score_type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "perlim-score-type", 
    "The type of scoring function to use. sp",
    (void *) &perlim_score_type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "mass-window", 
    "The peptide mass tolerance window", 
    (void *) &mass_window, 
    DOUBLE_ARG);
  
  /* Define required command line arguments */
  parse_arguments_set_req(
    "ms2", 
    "The name of the file (in MS2 format) from which to parse the spectrum.",
    (void *) &ms2_file, 
    STRING_ARG);

  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the MS2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within a .ms2 file.",
    (void *) &scan_num, 
    INT_ARG);

  parse_arguments_set_req(
    "fasta-file", 
    "The name of the file (in fasta format) from which to retrieve proteins and peptides.",
    (void *) &fasta_file, 
    STRING_ARG);
  
  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    //parse arguments
    SCORER_TYPE_T main_score = SP; 
    SCORER_TYPE_T perlim_score = SP; 
    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
    MATCH_COLLECTION_T* match_collection = NULL;
    MATCH_ITERATOR_T* match_iterator = NULL;
    MATCH_T* match = NULL;
    unsigned int max_rank = 500;

    //set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity", "verbosity level must be between 0-100");
    }

    //set verbosity
    set_verbosity_level(verbosity);

    //parse and update parameters
    parse_update_parameters(parameter_file);
    
    //check charge 
    if( charge < 1 || charge > 3){
      wrong_command(NULL, "The peptide charge. 1|2|3"); ///FIXME
    }

    //score type
    if(strcmp(get_string_parameter_pointer("score-type"), "sp")== 0){
      main_score = SP;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "xcorr")== 0){
      main_score = XCORR;
    }
    else{
      wrong_command(score_type, "The type of scoring function to use. sp");
    }

    //score type
    if(strcmp(get_string_parameter_pointer("perlim-score-type"), "sp")== 0){
      perlim_score = SP;
    }
    else{
      wrong_command(perlim_score_type, "The type of perliminary scoring function to use. sp");
    }
    
    //set max number of matches to return
    if(main_score == SP){
      //return all matches
      max_rank = _MAX_NUMBER_PEPTIDES;
    }
    else if(main_score == XCORR){
      // keep top 500 from SP to xcorr
      max_rank = 500;
    }
    else{
      //default 500
      max_rank = 500;
    }

    //parameters are now confirmed, can't be changed
    parameters_confirmed();
 
    //print header 1/2
    fprintf(stdout, "# SPECTRUM FILE: %s\n", ms2_file);
    fprintf(stdout, "# PROTEIN DATABASE: %s\n", fasta_file);
    fprintf(stdout, "# SPECTRUM SCAN NUMBER: %d\n", scan_num);


    //read ms2 file
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();

    //search for spectrum with correct scan number
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "no matching scan number: %d, in ms2 file: %s", scan_num, ms2_file);
      //free, exit
      exit(1);
    }

    //print header 2/2
    fprintf(stdout, "# SPECTRUM ID NUMBER: %d\n", get_spectrum_id(spectrum));
    fprintf(stdout, "# SPECTRUM PRECURSOR m/z: %.2f\n", get_spectrum_precursor_mz(spectrum));
    fprintf(stdout, "# SPECTRUM CHARGE: %d\n", charge);

    
    //get match collection with perlim match collection
    match_collection = new_match_collection_spectrum(spectrum, charge, max_rank, main_score);
    
    //later add scoring for main_score here!

    //create match iterator, TRUE: return match in sorted order of main_score type
    match_iterator = new_match_iterator(match_collection, main_score, TRUE);

    //print header
    if(main_score == SP){
      fprintf(stdout, "# %s\t%s\t%s\t%s\n", "sp_rank", "mass", "sp", "sequence");  
    }
    else if( main_score == XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "xcorr_rank", "sp_rank", "mass", "xcorr", "sp", "sequence");  
    }

    //iterate over all matches
    while(match_iterator_has_next(match_iterator)){
      match = match_iterator_next(match_iterator);
      print_match(match, stdout, TRUE, main_score);
    }

    //free match iterator
    free_match_iterator(match_iterator);
    free_match_collection(match_collection);
    free_spectrum_collection(collection);
    free_spectrum(spectrum);
  }
  else{
    char* usage = parse_arguments_get_usage("search_spectrum");
    result = parse_arguments_get_error(&error_message);
    fprintf(stderr, "Error in command line. Error # %d\n", result);
    fprintf(stderr, "%s\n", error_message);
    fprintf(stderr, "%s", usage);
    free(usage);
  }
  exit(0);
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
