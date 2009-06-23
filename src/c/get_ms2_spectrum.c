/*************************************************************************//**
 * \file get_ms2_spectrum.c
 * AUTHOR: Chris Park
 * CREATE DATE: 30 June 2006
 * DESCRIPTION: searches a given ms2 file for the spectrum with the given
 *              scan number.
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "parameter.h"
#include "carp.h"
#include "spectrum.h"
#include "peak.h"
#include "spectrum_collection.h"
#include "unistd.h"

#define NUM_MS2_OPTIONS 3
#define NUM_MS2_ARGUMENTS 3

int main(int argc, char** argv){

  /* Declarations */
  SPECTRUM_COLLECTION_T * collection; ///<spectrum collection
  SPECTRUM_T * spectrum; ///<the spectrum of interest
  BOOLEAN_T spectrum_found; ///<did we find the spectrum?
  int scan_number; ///< the query scan number
  FILE* output_file; ///< output file name
  //OUTPUT_WRITE_MODE_T file_flag;
  BOOLEAN_T options = FALSE; ///< do we want options?

  /* Define optional command line arguments */
  int num_options = NUM_MS2_OPTIONS;
  char* option_list[NUM_MS2_OPTIONS] = { 
    "version", 
    "stats", 
    "verbosity" }; //out-file

  int num_arguments = NUM_MS2_ARGUMENTS;
  char* argument_list[NUM_MS2_ARGUMENTS] = {
    "scan number",
    "ms2 file",
    "output file"
  };

  /* for debugging parameter handling */
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  //set_verbosity_level(CARP_DETAILED_DEBUG);
  set_verbosity_level(CARP_ERROR);

  /* set up parameters and default values */
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments );

  /* Parse the command line, for now param file not available */
  parse_cmd_line_into_params_hash(argc, argv, "crux-get-ms2-spectrum");

  /* Set verbosity */
  //set_verbosity_level(get_int_parameter("verbosity"));

  /* Get arguments */
  scan_number = get_int_parameter("scan number");
  char* ms2_filename = get_string_parameter_pointer("ms2 file");
  char* output_filename = get_string_parameter_pointer("output file");
  carp(CARP_DETAILED_DEBUG, "ms2_filename: %s", ms2_filename);
  carp(CARP_DETAILED_DEBUG, "output_filename: %s", output_filename);

  /* Get options */
  options = get_boolean_parameter("stats");
  //file_flag = get_output_mode_parameter("out-file");

  /* read input file */
  if( access(ms2_filename, F_OK)){
    carp(CARP_FATAL, "Could not read from ms2 file '%s'", ms2_filename);
  }
  carp(CARP_DETAILED_DEBUG, "Creating spectrum collection.");
  collection = new_spectrum_collection(ms2_filename);
  spectrum = allocate_spectrum();
  
  /* search for spectrum with correct scan number */
  spectrum_found = get_spectrum_collection_spectrum(collection, 
                                                    scan_number, 
                                                    spectrum);
  if( !spectrum_found ){
    carp(CARP_FATAL, "Could not find scan number %i", scan_number);
  }

  carp(CARP_INFO, "Found a spectrum with scan number %i", scan_number);

  if( !access(output_filename, F_OK)){//&& file_flag != FILE_REPLACE
    carp(CARP_FATAL, "The output file '%s' already exists.  " \
         "Use the --out-file option to overwrite or replace.",
         output_filename);
  }

  output_file = fopen(output_filename, "w");//output_mode_to_string(file_flag)

  if( access(output_filename, W_OK)){
    carp(CARP_FATAL, 
         "Could not write spectrum to '%s'.  Check file permissions", 
         output_filename);
  }

  print_spectrum(spectrum, output_file);
  fclose(output_file);

  int charge_state_index = 0; 
  int charge_state_num = get_spectrum_num_possible_z(spectrum);
  int* possible_z_array = get_spectrum_possible_z(spectrum);
  int possible_z;
  
  /* Print stats if requested */
  if(options){
    printf("Scan number: %i\n", scan_number);
    printf("Precursor m/z:%.2f\n", get_spectrum_precursor_mz(spectrum));
    printf("Total Ion Current:%.2f\n", get_spectrum_total_energy(spectrum));
    printf("Base Peak Intensity:%.1f\n", 
           get_spectrum_max_peak_intensity(spectrum)); // base is max
    printf("Number of peaks:%d\n", get_spectrum_num_peaks(spectrum));
    printf("Minimum m/z:%.1f\n", get_spectrum_min_peak_mz(spectrum));
    printf("Maximum m/z:%.1f\n", get_spectrum_max_peak_mz(spectrum));
    
    for(charge_state_index=0; charge_state_index < charge_state_num; 
                                                   ++charge_state_index){
      possible_z = possible_z_array[charge_state_index];
      printf("Charge state:%d\n", possible_z);
      printf("Neutral mass:%.2f\n", 
             get_spectrum_neutral_mass(spectrum, possible_z));
      printf("Charged mass:%.2f\n", get_spectrum_mass(spectrum, possible_z));
      printf("M+H+ mass:%.2f\n", 
             get_spectrum_singly_charged_mass(spectrum, possible_z));
    }
  }

  free(possible_z_array);
  free_spectrum(spectrum);
  free_spectrum_collection(collection);
  
  carp(CARP_INFO, "crux-get-ms2-spectrum finished.");
  return(0);
}
