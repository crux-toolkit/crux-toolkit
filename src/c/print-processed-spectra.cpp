/**
 * \file print-processed-spectra.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: September 18, 2009
 * \brief Main method for the print-processed-spectra command.
 *
 * For every spectrum in an ms2 file, process as for xcorr and print
 * peaks in ms2 format to new file.
 */

#include "print-processed-spectra.h"

int print_processed_spectra_main(int argc, char** argv){

  // Define optional command line arguments
  const char* option_list[] = { 
    "verbosity",
    "parameter-file", 
    "overwrite"
  };
  int num_options = sizeof(option_list) / sizeof(char*);

  // Define required command line arguments
  const char* argument_list[] = { "ms2 file", 
                                              "output file"}; 
  int num_arguments = sizeof(argument_list) / sizeof(char*);

  initialize_run(PROCESS_SPEC_COMMAND, 
                 argument_list, 
                 num_arguments,
                 option_list, num_options,
                 argc, argv);

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
  SpectrumCollection* spectra = new SpectrumCollection(input_ms2_name);
  if( spectra == NULL ){
    carp(CARP_FATAL, "Could not read spectra from %s.", input_ms2_name);
  }

  spectra->parse();
  carp(CARP_DEBUG, "Found %d spectra in file.", 
       spectra->getNumSpectra());

  // write header to output file
  char* header = spectra->getComment();
  fprintf(output_ms2, "%s", header);
  fprintf(output_ms2, "H\tComment\tSpectra processed as for Xcorr\n");

  // create iterator for getting spectra
  FilteredSpectrumChargeIterator* spectrum_iterator =
    new FilteredSpectrumChargeIterator(spectra);

  if( spectrum_iterator == NULL ){
    carp(CARP_FATAL, "Could create spectrum iterator");
  }

  // loop over all spectra, process, print
  while(spectrum_iterator->hasNext()){
    SpectrumZState cur_zstate;
    int cur_charge = 0;
    Spectrum* cur_spectrum = 
      spectrum_iterator->next(cur_zstate);

    cur_charge = cur_zstate.getCharge();
    carp(CARP_DETAILED_INFO, "Processing spectrum %d charge %d.",
         cur_spectrum->getFirstScan(), cur_charge);

    // change the peak values
    FLOAT_T* intensities = NULL;
    int max_mz_bin = 0;
    get_processed_peaks(cur_spectrum, cur_charge, XCORR,
                        &intensities, &max_mz_bin);

    // print processed spectrum
    cur_spectrum->printProcessedPeaks(cur_zstate, 
                                        intensities, max_mz_bin,
                                        output_ms2);
  }

  // close output file
  delete spectra;
  fclose(output_ms2);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux print-processed-spectra.");

  return(0);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

