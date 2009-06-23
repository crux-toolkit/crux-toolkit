/*************************************************************************//**
 * \file score_peptide_spectrum
 * AUTHOR: Chris Park
 * CREATE DATE: 10/13 2006
 * DESCRIPTION: Object for given a peptide and a spectrum, generate a
 * preliminary score(ex, Sp) 
 *
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
#include "parse_arguments.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "ion.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "scorer.h"

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){

  // print comment if given
  if(comment == NULL){
    carp(CARP_FATAL, "incorrect argument: %s", arg);
  }
  else {
    carp(CARP_FATAL, "incorrect argument: %s\n%s", arg, comment);
  }

}

int main(int argc, char** argv){

  // required variables
  char* ms2_file = NULL;
  int scan_num = 0;
  char* peptide_sequence = NULL;

  // optional variables
  char* charge = "2";
  char* type = "xcorr";
  char* parameter_file = "crux.params";
  int  verbosity = CARP_ERROR;

  // parsing variables
  int result = 0;
  char * error_message; 


  /* Define optional command line arguments */ 
  parse_arguments_set_opt(
    "charge", 
    "The peptide charge. 1|2|3",
    (void *) &charge, 
    STRING_ARG);
  
  parse_arguments_set_opt(
    "score-type", 
    "The type of scoring function to use. sp | xcorr",
    (void *) &type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "parameter-file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);

  parse_arguments_set_opt(
    "verbosity", 
    "Specify the verbosity of the current processes from 0-100.",
    (void *) &verbosity, 
    INT_ARG);

 /* Define required command line arguments */
  parse_arguments_set_req(
    "peptide-sequence", 
    "The literal peptide sequence (e.g. EAMAPK) that is used to predict the ions.", 
    (void *) &peptide_sequence, 
    STRING_ARG);

  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {
    // parsed arguments
    int peptide_charge = 1;
    SCORER_TYPE_T score_type = XCORR; 
    
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    ION_SERIES_T* ion_series = NULL;
    SCORER_T* scorer = NULL;
    FLOAT_T score = 0;
    
    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity", "verbosity level must be between 0-100");
    }

    // set verbosity
    set_verbosity_level(verbosity);

    // parse and update parameters
    parse_update_parameters(parameter_file);
    
    
    peptide_charge = get_int_parameter("charge");
    
    if( peptide_charge < 1 || peptide_charge > 3){
      wrong_command(charge, "The peptide charge. 1|2|3");
    }

    // check peptide sequence
    if(!valid_peptide_sequence(peptide_sequence)){
      wrong_command(peptide_sequence, "not a valid peptide sequence");
    }
    
    // score type
    if(strcmp(get_string_parameter_pointer("score-type"), "sp")== 0){
      score_type = SP;
    }
    else if(strcmp(get_string_parameter_pointer("score-type"), "xcorr")== 0){
      score_type = XCORR;
    }
    else{
      wrong_command(type, "The type of scoring function to use. sp | xcorr");
    }
    
    // parameters are now confirmed, can't be changed
    parameters_confirmed();
    
    // set ion constraint to sequest settings
    ION_CONSTRAINT_T* ion_constraint = NULL;
    
    if(score_type == SP){
      ion_constraint = new_ion_constraint_sequest_sp(peptide_charge);  
      // create new scorer
      scorer = new_scorer(SP);  
    }
    else if(score_type == XCORR){
      ion_constraint = new_ion_constraint_sequest_xcorr(peptide_charge);  
      scorer = new_scorer(XCORR);  
    }

    // create new ion series
    ion_series = new_ion_series(peptide_sequence, peptide_charge, ion_constraint);
   
   // now predict ions
    predict_ions(ion_series);
   
   // read ms2 file
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    
    // search for spectrum with correct scan number
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_FATAL, "failed to find spectrum with  scan_num: %d", scan_num);
    }
        
    // calculates the Sp score
    score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
   
    // print the Sp score
    if(score_type == SP){
      printf("Sp score is: %.2f\n", score);
    }
    else if(score_type == XCORR){
      printf("Xcorr score is: %.2f\n", score);
    }
    else{
      carp(CARP_ERROR, "invalid score type for the application");
    }
      
   // free heap
   free_scorer(scorer);
   free_ion_constraint(ion_constraint);
   free_ion_series(ion_series);
   free_spectrum_collection(collection);
   free_spectrum(spectrum);
   free_parameters();
 }
 else{
   char* usage = parse_arguments_get_usage("score_peptide_spectrum");
   result = parse_arguments_get_error(&error_message);
   carp(
     CARP_FATAL,
     "Error in command line. Error # %d\n%s\n%s",
     result,
     error_message,
     usage
   );
 }
 exit(0);
}
