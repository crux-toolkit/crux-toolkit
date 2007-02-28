/*****************************************************************************
 * \file search_spectra
 * AUTHOR: Chris Park
 * CREATE DATE: Jan 03 2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an optional parameter file, 
 * search all the spectrum against the peptides in the sequence database, and return high scoring peptides. 
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
#include "generate_peptides_iterator.h"
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
  char* usage = parse_arguments_get_usage("search_spectra");
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
  char* perlim_score_type = "sp";
  char* score_type = "sp";
  char* parameter_file = NULL;
  int verbosity = CARP_ERROR;
  
  //required
  char* ms2_file = NULL;
  char* fasta_file = NULL;
  double mass_window = 3;

  //parsing variables
  int result = 0;
  char * error_message;

  /* Define optional command line arguments */   
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
    "The type of scoring function to use. sp",
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
    SPECTRUM_ITERATOR_T* spectrum_iterator = NULL;
    MATCH_COLLECTION_T* match_collection = NULL;
    MATCH_ITERATOR_T* match_iterator = NULL;
    MATCH_T* match = NULL;
    int possible_charge = 0;
    int* possible_charge_array = NULL;
    int charge_index = 0;
    
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
    
    //score type
    if(strcmp(get_string_parameter_pointer("score-type"), "sp")== 0){
      main_score = SP;
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

    //always use index when search spectra!
    set_string_parameter("use-index", "T");

    //parameters are now confirmed, can't be changed
    parameters_confirmed();
 
    //print header
    fprintf(stdout, "# SPECTRUM FILE: %s\n", ms2_file);
    fprintf(stdout, "# PROTEIN DATABASE: %s\n", fasta_file);
    
    //read ms2 file
    collection = new_spectrum_collection(ms2_file);
    
    //parse the ms2 file for spectra
    if(!parse_spectrum_collection(collection)){
      carp(CARP_ERROR, "failed to parse ms2 file: %s", ms2_file);
      //free, exit
      exit(1);
    }
    
    //create spectrum iterator
    spectrum_iterator = new_spectrum_iterator(collection);
    
    //create a generate peptide iterator
    GENERATE_PEPTIDES_ITERATOR_T* peptide_iterator =  //NULL;//FIXME use neutral_mass, might chage to pick
      new_generate_peptides_iterator_mutable();


    //iterate over all spectrum in ms2 file
    while(spectrum_iterator_has_next(spectrum_iterator)){
      //get next spectrum
      spectrum = spectrum_iterator_next(spectrum_iterator);

      //get possible charge state
      possible_charge = get_spectrum_num_possible_z(spectrum);
      possible_charge_array = get_spectrum_possible_z_pointer(spectrum);
      
      //print spectrum info
      fprintf(stdout, "# SPECTRUM SCAN NUMBER: %d\n", get_spectrum_first_scan(spectrum));
      fprintf(stdout, "# SPECTRUM ID NUMBER: %d\n", get_spectrum_id(spectrum));
      fprintf(stdout, "# SPECTRUM PRECURSOR m/z: %.2f\n", get_spectrum_precursor_mz(spectrum));

      //iterate over all possible charge states
      for(charge_index = 0; charge_index < possible_charge; ++charge_index){

        //print working state
        fprintf(stdout, "# SPECTRUM CHARGE: %d\n", possible_charge_array[charge_index]);
        fprintf(stdout, "# %s\t%s\t%s\t%s\n", "rank", "mass", "sp", "sequence");  

        //DEBUG                
        //peptide_iterator =  //FIXME use neutral_mass, might chage to pick
        //new_generate_peptides_iterator_mutable();


        //get match collection with perlim match collection
        match_collection =
          new_match_collection_spectrum_with_peptide_iterator(spectrum, 
                                                              possible_charge_array[charge_index], 
                                                              500, perlim_score, peptide_iterator);
        
        
        //later add scoring for main_score here!
        
        //create match iterator, TRUE: return match in sorted order of main_score type
        match_iterator = new_match_iterator(match_collection, main_score, TRUE);
        
        //iterate over all matches
        while(match_iterator_has_next(match_iterator)){
          match = match_iterator_next(match_iterator);
          print_match(match, stdout, TRUE, SP);
        }
                
        //free match iterator
        free_match_iterator(match_iterator);
        free_match_collection(match_collection);

      }

      //THEy are all freed in free_spectrum collection
      //done with all charge states of this spectrum, Free!
      //free_spectrum(spectrum);
    }

    free_spectrum_iterator(spectrum_iterator);
    free_spectrum_collection(collection);
    //FIXME!!!
    //free up peptide iterator
    free_generate_peptides_iterator(peptide_iterator);
  }
  else{
    char* usage = parse_arguments_get_usage("search_spectra");
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
