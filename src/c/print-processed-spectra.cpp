/**
 * \file print-processed-spectra.cpp
 *
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 * REVISION:
 */

#include "print-processed-spectra.h"
#define NUM_PPS_OPTIONS 4
#define NUM_PPS_ARGS 2

int print_processed_spectra_main(int argc, char** argv){

  // Define optional command line arguments
  int num_options = NUM_PPS_OPTIONS;
  const char* option_list[NUM_PPS_OPTIONS] = { 
    "version",
    "verbosity",
    "parameter-file", 
    "overwrite"
  };

  // Define required command line arguments
  int num_arguments = NUM_PPS_ARGS ;
  const char* argument_list[NUM_PPS_ARGS] = { "ms2 file", 
                                              "output file"}; 

  // For output of parameter parsing
  set_verbosity_level(CARP_ERROR);  

  // set up parameters and their defaults in parameter.c
  initialize_parameters();

  // Define optional and required command line arguments
  select_cmd_line_options( option_list, num_options );
  select_cmd_line_arguments( argument_list, num_arguments);

  // Parse the command line, including the optional params file
  parse_cmd_line_into_params_hash(argc, argv, "crux print-processed-spectra");

  // Get arguments and options
  const char* input_ms2_name  = get_string_parameter_pointer("ms2 file");
  char* output_ms2_name = get_string_parameter("output file");
  prefix_fileroot_to_name(&output_ms2_name);
  const char* output_dir = get_string_parameter_pointer("output-dir");
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");

  // open output file
  create_output_directory(output_dir, overwrite);
  FILE* output_ms2 = create_file_in_path(output_ms2_name,
                                         output_dir,
                                         overwrite);
  // open input file
  SPECTRUM_COLLECTION_T* spectra = new_spectrum_collection(input_ms2_name);
  if( spectra == NULL ){
    carp(CARP_FATAL, "Could not read spectra from %s.", input_ms2_name);
  }

  parse_spectrum_collection(spectra);
  carp(CARP_DEBUG, "Found %d spectra in file.", 
       get_spectrum_collection_num_spectra(spectra));

  // write header to output file
  char* header = get_spectrum_collection_comment(spectra);
  fprintf(output_ms2, header);
  fprintf(output_ms2, "H\tComment\tSpectra processed as for Xcorr\n");

  // create iterator for getting spectra
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* spectrum_iterator =
    new_filtered_spectrum_charge_iterator(spectra);

  if( spectrum_iterator == NULL ){
    carp(CARP_FATAL, "Could create spectrum iterator");
  }

  // loop over all spectra, process, print
  while(filtered_spectrum_charge_iterator_has_next(spectrum_iterator)){
    int cur_charge = 0;
    SPECTRUM_T* cur_spectrum = 
      filtered_spectrum_charge_iterator_next(spectrum_iterator, &cur_charge);

    carp(CARP_DETAILED_INFO, "Processing spectrum %d charge %d.",
         get_spectrum_first_scan(cur_spectrum), cur_charge);

    // change the peak values
    FLOAT_T* intensities = NULL;
    int max_mz_bin = 0;
    get_processed_peaks(cur_spectrum, cur_charge,
                        &intensities, &max_mz_bin);

    // print processed spectrum
    print_spectrum_processed_peaks(cur_spectrum, cur_charge, 
                                   intensities, max_mz_bin,
                                   output_ms2);
  }

  // close output file
  free_spectrum_collection(spectra);
  fclose(output_ms2);

  carp(CARP_INFO, "Finished processing spectra.");
  
  return 0;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

