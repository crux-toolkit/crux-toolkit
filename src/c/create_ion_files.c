/*****************************************************************************
 * \file create_ion_files
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 8/8 2007
 * DESCRIPTION: Creates files describing ion series, for input to GMTK.
 * REVISION: $Revision: 1.4 $
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
  char* usage = parse_arguments_get_usage("predict_peptide_ions");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  //print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

int main(int argc, char** argv){

  //required variables
  char* ms2_file = NULL;
  int scan_num = 0;
  char* peptide_sequence = NULL;
  char* output_directory = NULL;

  //optional variables
  char* charge = "2";
  char* parameter_file = NULL;
  int  verbosity = CARP_ERROR;

  //parsing variables
  int result = 0;
  char * error_message; 


  /* Define optional command line arguments */ 
  parse_arguments_set_opt(
    "charge", 
    "The peptide charge. 1|2|3",
    (void *) &charge, 
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
    "The peptide sequence (e.g. EAMAPK) that is used to predict the ions.", 
    (void *) &peptide_sequence, 
    STRING_ARG);

  parse_arguments_set_req(
    "scan-number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2-filename", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);

	parse_arguments_set_req(
    "output-dir", 
    "A directory in which to place the ion files.",
    (void *) &output_directory,
    STRING_ARG);

  /* Parse the command line */
  if (parse_arguments(argc, argv, 0)) {

    //parsed arguments
    int peptide_charge = 2;

    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    ION_SERIES_T* ion_series = NULL;
    
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
    
    peptide_charge = get_int_parameter("charge");
    
    if(peptide_charge < 1 || peptide_charge > 3){
      wrong_command(charge, "The peptide charge. 1|2|3");
    }

    //check peptide sequence
    if(!valid_peptide_sequence(peptide_sequence)){
      wrong_command(peptide_sequence, "not a valid peptide sequence");
    }
    
    //parameters are now confirmed, can't be changed
    parameters_confirmed();
    
    //read ms2 file
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    
    //search for spectrum with correct scan number
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "failed to find spectrum with  scan_num: %d", scan_num);
      free_spectrum_collection(collection);
      free_spectrum(spectrum);
      exit(1);
    }
    
    //set ion constraint to sequest settings
    ION_CONSTRAINT_T* ion_constraint = new_ion_constraint_gmtk(peptide_charge);  
    //create new ion series
    ion_series = new_ion_series(
									 peptide_sequence, peptide_charge, ion_constraint);
   
   	//now predict ions
   	predict_ions(ion_series);
       
	 	// TODO figure out how to handle neutral losss
	 	// TODO figure out the GMTK file input lists
   	// GMTK output peptide ion files
   	if (output_ion_files(output_directory, spectrum, ion_series) == FALSE){
  		 carp(CARP_FATAL, "failed to create ion files for: %s %i %s", 
				 ms2_file, scan_num, peptide_sequence);
	 	}
   
   	//free heap
   	free_ion_constraint(ion_constraint);
   	free_ion_series(ion_series);
   	free_spectrum_collection(collection);
   	free_spectrum(spectrum);
   	free_parameters();
 	}
 	else{
   	char* usage = parse_arguments_get_usage("create_ion_files");
   	result = parse_arguments_get_error(&error_message);
   	fprintf(stderr, "Error in command line. Error # %d\n", result);
   	fprintf(stderr, "%s\n", error_message);
   	fprintf(stderr, "%s", usage);
   	free(usage);
 	}
 	exit(0);
}
