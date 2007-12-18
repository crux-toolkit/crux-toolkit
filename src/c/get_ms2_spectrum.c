/*****************************************************************************
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

#define NUM_MS2_OPTIONS 2
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
  /*
  char* USAGE = 
    "USAGE: get_ms2_spectrum"
    "\t<scan number>\n"
    "\t<input filename>\n"
    "\t<output filename>\n"
    "\t[--stats]\n"
    "SYNOPSIS: get_ms2_spectrum scan_number input_file input_file [--stats]\n"  
    "OPTIONS: --stats\n"
    "\tOutput some summary statistics for the particular spectrum to STDOUT\n"; 
  */
  /* Define optional command line arguments */
  int num_options = NUM_MS2_OPTIONS;
  char* option_list[NUM_MS2_OPTIONS] = { "stats", "verbosity" }; //out-file

  int num_arguments = NUM_MS2_ARGUMENTS;
  char* argument_list[NUM_MS2_ARGUMENTS] = {
    "scan number",
    "ms2 file",
    "output file"
  };

  /* for debugging parameter handling */
  //set_verbosity_level(CARP_DETAILED_DEBUG);

  /* set up parameters and default values */
  initialize_parameters();

  /* Define optional and required arguments */
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments );

  /* Parse the command line, for now param file not available */
  parse_cmd_line_into_params_hash(argc, argv);

  /* Set verbosity */
  set_verbosity_level(get_int_parameter("verbosity"));

  /*
  // check command line argument count
  if (argc != 4) {
    if(argc == 5 && 
       strcmp(argv[4], "--stats")==0){
      options = TRUE;      
    }
    else{
      fprintf (stderr, "%s", USAGE);
      exit(1);
    }
  }
  
  // check if first argument is type int for scan number
  if (sscanf(argv[1], "%d", &scan_num) != 1) {
    fprintf (stderr, "ERROR:second arguement must be type int\n");
    fprintf (stderr, "%s", USAGE);
    exit(1);
  }

  // check if the input file exist
  if(access(argv[2], F_OK)){
    fprintf(stderr,"input file:\"%s\" could not be opened\n", argv[2]);
    exit(1);
  }
  
  // check if the output file doesn't aleady exist
  if(!access(argv[3], F_OK)){
    fprintf(stderr,"output file:\"%s\" already exist\n", argv[3]);
    exit(1);
  }
  */

  /* Get arguments */
  scan_number = get_int_parameter("scan number");
  char* ms2_filename = get_string_parameter_pointer("ms2 file");
  char* output_filename = get_string_parameter_pointer("output file");

  /* Get options */
  options = get_boolean_parameter("stats");
  //file_flag = get_output_mode_parameter("out-file");


  /* read input file */
  //  collection = new_spectrum_collection(argv[2]);
  if( access(ms2_filename, F_OK)){
    carp(CARP_FATAL, "Could not read from ms2 file '%s'", ms2_filename);
    exit(1);
  }
  collection = new_spectrum_collection(ms2_filename);
  spectrum = allocate_spectrum();
  
  /* search for spectrum with correct scan number */
  spectrum_found = get_spectrum_collection_spectrum(collection, 
							  scan_number, 
							  spectrum);
//  if(get_spectrum_collection_spectrum(collection, scan_number, spectrum)){ 
  if( !spectrum_found ){
    carp(CARP_FATAL, "Could not find scan number %i", scan_number);
    exit(1);
  }

  //    printf("Found a spectrum with correct scan number\n");
  carp(CARP_INFO, "Found a spectrum with scan number %i", scan_number);
  //    output_file = fopen(argv[3], "w");

  if( !access(output_filename, F_OK)){//&& file_flag != FILE_REPLACE
    carp(CARP_FATAL, "The output file '%s' already exists.  " \
	 "Use the --out-file option to overwrite or replace.",
	 output_filename);
    exit(1);
  }

  output_file = fopen(output_filename, "w");//output_mode_to_string(file_flag)

  if( access(output_filename, W_OK)){
    carp(CARP_FATAL, 
	 "Could not write spectrum to '%s'.  Check file permissions", 
	 output_filename);
    exit(1);
  }

  print_spectrum(spectrum, output_file);
  fclose(output_file);
  //spectrum_found = TRUE;

/*  else{
    fprintf(stderr,
            "Warning:\n" 
            "The ms2 file %s does not contain the spectrum with scan number %d.\n", 
            argv[2], 
            scan_number);
    spectrum_found = FALSE;
  }
*/
  int charge_state_index = 0; 
  int charge_state_num = get_spectrum_num_possible_z(spectrum);
  int* possible_z_array = get_spectrum_possible_z(spectrum);
  int possible_z;
  
  // print only if found spectrum and with option flag
  /* Print stats if requested */
  //  if(spectrum_found && options){
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
  
  // failed to find spectrum
  /*
  if(!spectrum_found){
    return(1);
  }
  */
  // exit success

  carp(CARP_INFO, "crux-get-ms2-spectrum finished.");
  return(0);
}
