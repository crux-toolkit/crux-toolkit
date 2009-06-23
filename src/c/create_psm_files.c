/*************************************************************************//**
 * \file create_ion_files
 * AUTHOR: Aaron Klammer
 * CREATE DATE: 8/8 2007
 * DESCRIPTION: Creates files describing ion series, for input to GMTK.
 * REVISION: $Revision: 1.11 $
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

  // print comment if given
  if(comment == NULL){
    carp(CARP_FATAL, "incorrect argument: %s", arg);
  }
  else {
    carp(CARP_FATAL, "incorrect argument: %s\n%s", arg, comment);
  }
}

#define PSM_NUM_OPTIONS 4
#define PSM_NUM_ARGUMENTS 5

int main(int argc, char** argv){

  // required variables
  char* ms2_file;
  int scan_num;
  char* peptide_file_name;
  char* output_directory;

  // optional variables
  int charge;
  int starting_sentence_idx;
  char* model_type;

  /* Define optional command line arguments */ 
  int num_options = PSM_NUM_OPTIONS;
  char* option_list[PSM_NUM_OPTIONS] = { 
    "verbosity",
    "parameter-file", 
    "charge",
    "starting-sentence-idx"
  };

  /* Define required command line arguments */ 
  int num_arguments = PSM_NUM_ARGUMENTS;
  char* argument_list[PSM_NUM_ARGUMENTS] = { 
    "peptide-file-name",
    "scan number",
    "ms2 file",
    "output-dir",
    "model-type"
  };

  set_verbosity_level(CARP_ERROR);
  carp(CARP_DETAILED_DEBUG, "Creating ion files");


  initialize_parameters();
  select_cmd_line_options(option_list, num_options);
  select_cmd_line_arguments(argument_list, num_arguments);

  parse_cmd_line_into_params_hash(argc, argv, "create_psm_files");

  // parsed arguments
  ms2_file = get_string_parameter("ms2 file");
  scan_num = get_int_parameter("scan number");
  peptide_file_name = get_string_parameter("peptide-file-name");
  output_directory = get_string_parameter("output-dir");
  model_type = get_string_parameter("model-type");
  charge = get_int_parameter("charge");
  starting_sentence_idx = get_int_parameter("starting-sentence-idx");

  SPECTRUM_T* spectrum = NULL;
  SPECTRUM_COLLECTION_T* collection = NULL;
 
  // read ms2 file
  carp(CARP_INFO, "Reading ms2 file %s", ms2_file);
  collection = new_spectrum_collection(ms2_file);
  spectrum = allocate_spectrum();
  
  // search for spectrum with correct scan number
  carp(CARP_INFO, "Retrieving spectrum %i", scan_num);
  if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
    free_spectrum_collection(collection);
    free_spectrum(spectrum);
    carp(CARP_ERROR, "Failed to find spectrum with scan_num: %d", scan_num);
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
      
  carp(CARP_INFO, "Creating and outputting ions");

  // output GMTK peptide ion files
  // PAIRED add output_psm_files_paired and change this to 
  // output_psm_files_single
  BOOLEAN_T return_status = FALSE;
  if (strcmp(model_type, "single") == 0){
    return_status = output_psm_files_single(
        output_directory, spectrum, peptides, num_lines, 
        charge, starting_sentence_idx);
  } else if (strcmp(model_type, "paired") == 0){
    return_status = output_psm_files_paired(
        output_directory, spectrum, peptides, num_lines, 
        charge, starting_sentence_idx);
  }

  if (return_status == FALSE){
    carp(CARP_FATAL, "Failed to create ion files for: %s %i %s.", 
       ms2_file, scan_num, peptide_file_name);
  } else {
    carp(CARP_INFO, "Done outputting files.");
  }

  // free heap 
  int peptide_idx = 0;
  for (peptide_idx=0; peptide_idx < num_lines; peptide_idx++){
    free(peptides[peptide_idx]);
  }
  free(peptides);
  free_spectrum_collection(collection);
  free_spectrum(spectrum);
  free_parameters();
  exit(0);
}
