/*****************************************************************************
 * \file create_ion_files
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 8/8 2007
 * DESCRIPTION: Creates files describing ion series, for input to GMTK.
 * REVISION: $Revision: 1.17 $
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

#define MAX_PEPTIDES 5000

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
  char* peptide_file_name = NULL;
  char* output_directory = NULL;

  //optional variables
  int charge = 2;
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
    INT_ARG);
  
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
    "peptide-file-name", 
    "A file containing the peptide sequences for which to create ions files.", 
    (void *) &peptide_file_name, 
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

    // parsed arguments
    SPECTRUM_T* spectrum = NULL;
    SPECTRUM_COLLECTION_T* collection = NULL;
    ION_SERIES_T* ion_series = NULL;
    
    // set verbosity
    if(CARP_FATAL <= verbosity && verbosity <= CARP_MAX){
      set_verbosity_level(verbosity);
    }
    else{
      wrong_command("verbosity", "verbosity level must be between 0-100");
    }

    // parse and update parameters
    parse_update_parameters(parameter_file);
    
    // parameters now can't be changed
    parameters_confirmed();
	
    // read ms2 file
    carp(CARP_INFO, "Reading ms2 file %s", ms2_file);
    collection = new_spectrum_collection(ms2_file);
    spectrum = allocate_spectrum();
    
    // search for spectrum with correct scan number
    carp(CARP_INFO, "Retrieving spectrum %i", scan_num);
    if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
      carp(CARP_ERROR, "Failed to find spectrum with scan_num: %d", scan_num);
      free_spectrum_collection(collection);
      free_spectrum(spectrum);
      exit(1);
    }

		// prepare the spectrum 
    carp(CARP_INFO, "Normalizing spectrum %i", scan_num);
		sum_normalize_spectrum(spectrum);

    carp(CARP_INFO, "Ranking spectrum peaks %i", scan_num);
    spectrum_rank_peaks(spectrum); 

		// parse the peptides
		int num_lines;
    carp(CARP_INFO, "Parsing peptides from %s", peptide_file_name);
		char** peptides = parse_file(peptide_file_name, MAX_PEPTIDES, &num_lines);
				
		int peptide_idx = 0;
		char* peptide_sequence = NULL;
    carp(CARP_INFO, "Creating and outputting ions");

    // output GMTK peptide ion files
    if (output_psm_files(output_directory, spectrum, peptides, num_lines) 
          == FALSE){
       carp(CARP_FATAL, "Failed to create ion files for: %s %i %s.", 
          ms2_file, scan_num, peptide_file_name);
    } else {
      carp(CARP_INFO, "Done outputting files.");
    }

    // free heap 
    int peptide_idx;
    for (peptide_idx=0; peptide_idx < num_lines; peptide_idx++){
      free(peptides[peptide_idx]);
    }
    free(peptides);
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
