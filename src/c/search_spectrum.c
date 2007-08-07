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
  char* prelim_score_type = "sp";
  char* score_type = "xcorr";
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
    "The type of scoring function to use: "
    "logp_exp_sp | logp_bonf_exp_sp | logp_evd_xcorr | logp_bonf_evd_xcorr |"
    "logp_weibull_sp | logp_bonf_weibull_sp | logp_weibull_xcorr | "
    "logp_bonf_weibull_xcorr | xcorr",
    (void *) &score_type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "prelim-score-type", 
    "The type of scoring function to use. sp",
    (void *) &prelim_score_type, 
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
    SCORER_TYPE_T main_score = LOGP_EXP_SP; 
    SCORER_TYPE_T prelim_score = SP; 
    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
    MATCH_COLLECTION_T* match_collection = NULL;
    MATCH_ITERATOR_T* match_iterator = NULL;
    MATCH_T* match = NULL;
    long int max_rank_preliminary = 500;
    long int max_rank_result = 500;
    float mass_offset = 0;

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
    
    //always use index when search spectrum!
    set_string_parameter("use-index", "T");
    
    //parameters are now confirmed, can't be changed
    parameters_confirmed();

    /***** All parameters must be taken through get_*_parameter() method ******/

    //main score type
    //score type
    char* score_type = get_string_parameter_pointer("score-type");
    if(strcmp(score_type, "logp_exp_sp")== 0){
      main_score = LOGP_EXP_SP;
    }
    else if(strcmp(score_type, "logp_bonf_exp_sp")== 0){
      main_score = LOGP_BONF_EXP_SP;
    }    
    else if(strcmp(score_type, "logp_weibull_sp")== 0){
      main_score = LOGP_WEIBULL_SP;
    }
    else if(strcmp(score_type, "logp_bonf_weibull_sp")== 0){
      main_score = LOGP_BONF_WEIBULL_SP;
    }    
    else if(strcmp(score_type, "logp_weibull_xcorr")== 0){
      main_score = LOGP_WEIBULL_XCORR;
    }
    else if(strcmp(score_type, "logp_bonf_weibull_xcorr")== 0){
      main_score = LOGP_BONF_WEIBULL_XCORR;
    }    
    else if(strcmp(score_type, "xcorr")== 0){
      main_score = XCORR;
    }
    else if(strcmp(score_type, "logp_evd_xcorr")== 0){
      main_score = LOGP_EVD_XCORR;
    }
    else if(strcmp(score_type, "logp_bonf_evd_xcorr")== 0){
      main_score = LOGP_BONF_EVD_XCORR;
    }
    else{
      wrong_command(score_type, "The type of scoring function must be: "
          "logp_exp_sp | logp_bonf_exp_sp | "
          "logp_weibull_sp | logp_bonf_weibull_sp |"
          "logp_weibull_xcorr | logp_bonf_weibull_xcorr |"
          "logp_evd_xcorr | logp_bonf_evd_xcorr | xcorr");
    }

    //preliminary score type
    //score type
    if(strcmp(get_string_parameter_pointer("prelim-score-type"), "sp")== 0){
      prelim_score = SP;
    }
    else{
      wrong_command(prelim_score_type, "The type of preliminary scoring function to use. sp");
    }
    
    //get mass offset from precursor mass to search for candidate peptides
    mass_offset = get_double_parameter("mass-offset", 0);
    
    //set max number of preliminary scored peptides to use for final scoring
    max_rank_preliminary = get_int_parameter("max-rank-preliminary", 500);

    //set max number of final scoring matches to print as output
    max_rank_result = get_int_parameter("max-rank-result", 500);
    
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
    fprintf(stdout, "# MASS OFFSET: %.2f\n", mass_offset);

    //get match collection with prelim match collection
    match_collection = new_match_collection_spectrum(spectrum, charge, max_rank_preliminary, 
                                                     prelim_score, main_score, mass_offset, FALSE);
    
    //create match iterator, TRUE: return match in sorted order of main_score type
    match_iterator = new_match_iterator(match_collection, main_score, TRUE);

    //print additional header
    fprintf(stdout, "# PEPTIDES SEARCHED: %d\n", get_match_collection_experimental_size(match_collection));
    fprintf(stdout, "# PEPTIDES SAMPLED FOR LOGP_SP: %d\n",get_match_collection_top_fit_sp(match_collection) );
    
    //print header
    if(main_score == LOGP_EXP_SP){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_exp_sp_rank", "sp_rank", "mass", "logp_exp_sp", "sp", "sequence");  
    }
    else if(main_score == LOGP_BONF_EXP_SP){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_bonf_exp_sp_rank", "sp_rank", "mass", "logp_bonf_exp_sp", "sp", "sequence");  
    }
    else if(main_score == LOGP_BONF_WEIBULL_SP){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_bonf_weibull_sp_rank", "sp_rank", "mass", "logp_bonf_weibull_sp", "sp", "sequence");  
    }
    else if(main_score == LOGP_WEIBULL_SP){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_weibull_sp_rank", "sp_rank", "mass", "logp_weibull_sp", "sp", "sequence");  
    }
    else if(main_score == LOGP_BONF_WEIBULL_XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_bonf_weibull_xcorr_rank", "sp_rank", "mass", "logp_bonf_weibull_xcorr", "sp", "sequence");  
    }
    else if(main_score == LOGP_WEIBULL_XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_weibull_xcorr_rank", "sp_rank", "mass", "logp_weibull_xcorr", "sp", "sequence");  
    }
    else if(main_score == XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "xcorr_rank", "sp_rank", "mass", "xcorr", "sp", "sequence");  
    }
    else if(main_score == LOGP_EVD_XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_evd_xcorr_rank", "sp_rank", "mass", "logp_evd_xcorr", "sp", "sequence");  
    }
    else if(main_score == LOGP_BONF_EVD_XCORR){
      fprintf(stdout, "# %s\t%s\t%s\t%s\t%s\t%s\n", "logp_bonf_evd_xcorr_rank", "sp_rank", "mass", "logp_bonf_evd_xcorr", "sp", "sequence");  
    }
    

    //iterate over matches
    int match_count = 0;
    while(match_iterator_has_next(match_iterator)){
      ++match_count;
      match = match_iterator_next(match_iterator);
      print_match(match, stdout, TRUE, main_score);
      
      //print only up to max_rank_result of the matches
      if(match_count >= max_rank_result){
        break;
      }
    }

    //free match iterator
    free_match_iterator(match_iterator);
    free_match_collection(match_collection);
    free_spectrum_collection(collection);
    free_spectrum(spectrum);

    //test function
    /*
    pcEndRegistration();
    carp(CARP_DEBUG, "Percolator, did it happen?");
    */
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
