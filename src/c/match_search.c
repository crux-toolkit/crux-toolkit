/*****************************************************************************
 * \file match_search.c
 * AUTHOR: Chris Park
 * CREATE DATE: 6/18/2007
 * DESCRIPTION: Given as input an ms2 file, a sequence database, and an optional parameter file, 
 * search all the spectrum against the peptides in the sequence database, and return high scoring peptides. 
 * ouput as binary ouput and optional sqt file format
 * REVISION: 
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
  char* prelim_score_type = "sp";
  char* score_type = "xcorr";
  char* parameter_file = NULL;
  int verbosity = CARP_ERROR;
  char* match_output_folder = "."; //FIXME if needed..
  char* output_mode = "binary";
  char* sqt_output_file = "Prefix of <ms2 input filename>.psm";
  double spectrum_min_mass = 0;
  double spectrum_max_mass = INFINITY;
  int number_decoy_set = 2;
    
  //required
  char* ms2_file = NULL;
  char* fasta_file = NULL;
  
  //parsing variables
  int result = 0;
  char* error_message;

  /* Define optional command line arguments */   
  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);

  parse_arguments_set_opt(
    "match-output-folder",
    "The output folder to store the serialized binary output files.",
    (void *) &match_output_folder,
    STRING_ARG); 
  
  parse_arguments_set_opt(
    "output-mode",
    "The output format. sqt|binary|all",
    (void *) &output_mode,
    STRING_ARG); 
  
  parse_arguments_set_opt(
    "sqt-output-file",
    "Name of the sqt file to place matches.",
    (void *) &sqt_output_file,
    STRING_ARG); 
  
  parse_arguments_set_opt(
    "spectrum-min-mass", 
    "The lowest spectrum m/z to search in the ms2 file.",
    (void *) &spectrum_min_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "spectrum-max-mass", 
    "The highest spectrum m/z to search in the ms2 file.",
    (void *) &spectrum_max_mass, 
    DOUBLE_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG); 
  
  parse_arguments_set_opt(
    "score-type", 
    "The type of scoring function to use. logp_exp_sp | logp_bonf_exp_sp | logp_evd_xcorr | logp_bonf_evd_xcorr | xcorrlogp_exp_sp | logp_bonf_exp_sp | xcorr",
    (void *) &score_type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "prelim-score-type", 
    "The type of preliminary scoring function to use. sp",
    (void *) &prelim_score_type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "number-decoy-set", 
    "Specify the number of decoy sets to generate for match_analysis.",
    (void *) &number_decoy_set, 
    INT_ARG);
  
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
    SCORER_TYPE_T main_score = XCORR; 
    SCORER_TYPE_T prelim_score = SP; 
    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL; ///<spectrum collection
    SPECTRUM_ITERATOR_T* spectrum_iterator = NULL;
    MATCH_COLLECTION_T* match_collection = NULL;
    int possible_charge = 0;
    int* possible_charge_array = NULL;
    int charge_index = 0;
    long int max_rank_preliminary = 500;
    long int max_rank_result = 500;
    int top_match = 1;
    MATCH_SEARCH_OUPUT_MODE_T output_type = BINARY_OUTPUT;
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
    
    //always use index when search spectra!
    set_string_parameter("use-index", "T");

    //generate sqt ouput file if not set by user
    if(strcmp(get_string_parameter_pointer("sqt-output-file"), "Prefix of <ms2 input filename>.psm") ==0){
      sqt_output_file = generate_name(ms2_file, ".psm", ".ms2", NULL);
      set_string_parameter("sqt-output-file", sqt_output_file);
    }
    
    //parameters are now confirmed, can't be changed
    parameters_confirmed();
    
    /***** Now, must get all parameters through get_*_parameter ****/
    
    //number_decoy_set
    number_decoy_set = get_int_parameter("number-decoy-set", 2);

    //main score type
    if(strcmp(get_string_parameter_pointer("score-type"), "logp_exp_sp")== 0){
      main_score = LOGP_EXP_SP;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "logp_bonf_exp_sp")== 0){
      main_score = LOGP_BONF_EXP_SP;
    }    
    else if(strcmp(get_string_parameter_pointer("score-type"), "xcorr")== 0){
      main_score = XCORR;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "logp_evd_xcorr")== 0){
      main_score = LOGP_EVD_XCORR;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "logp_bonf_evd_xcorr")== 0){
      main_score = LOGP_BONF_EVD_XCORR;
    }
    else{
      wrong_command(score_type, "The type of scoring function to use. logp_exp_sp | logp_bonf_exp_sp | logp_evd_xcorr | logp_bonf_evd_xcorr | xcorr");
    }
    
    //preliminary score type
    if(strcmp(get_string_parameter_pointer("prelim-score-type"), "sp")== 0){
      prelim_score = SP;
    }
    else{
      wrong_command(prelim_score_type, "The type of preliminary scoring function to use. sp");
    }
        
    //get output-mode
    if(strcmp(get_string_parameter_pointer("output-mode"), "binary")== 0){
      output_type = BINARY_OUTPUT;
    }
    else if(strcmp(get_string_parameter_pointer("output-mode"), "sqt")== 0){
      output_type = SQT_OUTPUT;
    }
    else if(strcmp(get_string_parameter_pointer("output-mode"), "all")== 0){
      output_type = ALL_OUTPUT;
    }
    else{
      wrong_command(output_mode, "The output mode to use. binary|sqt|all");
    }
    
    //set max number of preliminary scored peptides to use for final scoring
    max_rank_preliminary = get_int_parameter("max-rank-preliminary", 500);

    //set max number of final scoring matches to print as output in sqt
    max_rank_result = get_int_parameter("max-rank-result", 500);

    //set max number of matches to be serialized per spectrum
    top_match = get_int_parameter("top-match", 1);

    //get mass offset from precursor mass to search for candidate peptides
    mass_offset = get_double_parameter("mass-offset", 0);    
    
    /************** done with parameter setting **************/
    
    char* psm_result_filename = NULL;
    FILE* psm_result_file = NULL;
    FILE* psm_result_file_sqt = NULL;
    
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
    
    //get psm_result file handler
    psm_result_file = 
      get_spectrum_collection_psm_result_filename(collection,
                                                  match_output_folder,
                                                  &psm_result_filename,
                                                  ".ms2"
                                                  );
    
    //get psm_result sqt file handle if needed
    if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
      psm_result_file_sqt = 
        create_file_in_path(sqt_output_file, match_output_folder);
    }
    
    //did we get the file handles?
    if(psm_result_file == NULL ||
       ((output_type == SQT_OUTPUT || output_type == ALL_OUTPUT) &&
       psm_result_file_sqt == NULL)){
      carp(CARP_ERROR, "failed with output mode");
      free(sqt_output_file);
      free(psm_result_filename);
      exit(-1);
    }
    
    //serialize the header information
    serialize_header(collection, fasta_file, psm_result_file);
    
    int spectra_idx = 0;
    //iterate over all spectrum in ms2 file
    while(spectrum_iterator_has_next(spectrum_iterator)){
      //get next spectrum
      spectrum = spectrum_iterator_next(spectrum_iterator);

      //select spectra that are within m/z target interval
      if(get_spectrum_precursor_mz(spectrum) < spectrum_min_mass ||
         get_spectrum_precursor_mz(spectrum) > spectrum_max_mass)
        {
          continue;
        }
      
      //get possible charge state
      possible_charge = get_spectrum_num_possible_z(spectrum);
      possible_charge_array = get_spectrum_possible_z_pointer(spectrum);
      
      //iterate over all possible charge states
      for(charge_index = 0; charge_index < possible_charge; ++charge_index){
        ++spectra_idx;
        
        //get match collection with scored, ranked match collection
        match_collection =
          new_match_collection_spectrum(spectrum, 
                                        possible_charge_array[charge_index], 
                                        max_rank_preliminary, prelim_score, 
                                        main_score, mass_offset, FALSE);
        
        //serialize the psm features to ouput file upto 'top_match' number of 
        //top peptides among the match_collection
        serialize_psm_features(match_collection, psm_result_file, top_match, prelim_score, main_score);
        
        //should I ouput the match_collection result as a SQT file?
        //FIXME ONLY one header
        if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
          print_match_collection_sqt(psm_result_file_sqt, max_rank_result,
                                     match_collection, spectrum, 
                                     prelim_score, main_score);
        }        
        free_match_collection(match_collection);
      }
    }
    
    free(sqt_output_file);
    if(output_type == SQT_OUTPUT || output_type == ALL_OUTPUT){
      fclose(psm_result_file_sqt);
    }
    fclose(psm_result_file);
    free(psm_result_filename);
    free_spectrum_iterator(spectrum_iterator);
    free_spectrum_collection(collection);
  }
  else{
    char* usage = parse_arguments_get_usage("match_search");
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
