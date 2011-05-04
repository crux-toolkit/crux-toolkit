/*************************************************************************//**
 * \file get_ms2_spectrum.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 30 June 2006
 * \brief searches a given ms2 file for the spectrum with the given
 * scan number.
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <vector>
#include "parameter.h"
#include "carp.h"
#include "Spectrum.h"
#include "Peak.h"
#include "SpectrumCollectionFactory.h"
#include "unistd.h"

/****************************************************************************
 * Read a string into either a single positive integer or a range of
 * positive integers, depending on whether it contains a hyphen.
 ****************************************************************************/
static void parse_scan_numbers
  (char* my_string,
   int*  min_scan,
   int*  max_scan)
{
  // Does the string contain a hyphen?
  int string_index = 0;
  while (my_string[string_index] != '\0') {
    if (my_string[string_index] == '-') {
      my_string[string_index] = '\0';
      *min_scan = atoi(my_string);
      *max_scan = atoi(my_string + string_index + 1);
      assert((*min_scan > 0) && (*max_scan > 0)); 
      assert(*min_scan < *max_scan);
      return;
    }
    string_index++;
  }

  // If no hyphen, then the range is just a single integer.
  *min_scan = *max_scan = atoi(my_string);
}

/****************************************************************************
 * MAIN
 ****************************************************************************/
static const int NUM_MS2_OPTIONS = 3;
static const int NUM_MS2_ARGUMENTS = 2;

int main(int argc, char** argv){

  /* Declarations */
  BOOLEAN_T options = FALSE; ///< do we want options?

  /* Define optional command line arguments */
  int num_options = NUM_MS2_OPTIONS;
  const char* option_list[NUM_MS2_OPTIONS] = { 
    "version", 
    "stats", 
    "verbosity"/*,
    "use-mstoolkit"*/};

  int num_arguments = NUM_MS2_ARGUMENTS;
  const char* argument_list[NUM_MS2_ARGUMENTS] = {
    "scan number",
    "ms2 file"
  };

  /* for debugging parameter handling */
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
  int min_scan;
  int max_scan;
  parse_scan_numbers(get_string_parameter("scan number"), &min_scan, &max_scan);
  fprintf(stderr, "Scanning from %d to %d.\n", min_scan, max_scan);

  const char* ms2_filename = get_string_parameter_pointer("ms2 file");
  carp(CARP_DETAILED_DEBUG, "ms2_filename: %s", ms2_filename);

  /* Get options */
  options = get_boolean_parameter("stats");

  /* read input file */
  if( access(ms2_filename, F_OK)){
    carp(CARP_FATAL, "Could not read from ms2 file '%s'", ms2_filename);
  }
  carp(CARP_DETAILED_DEBUG, "Creating spectrum collection.");
  SpectrumCollection* collection = SpectrumCollectionFactory::create(ms2_filename);

  int num_found = 0;
  for (int scan_number = min_scan; scan_number <= max_scan; scan_number++) {

    /* search for spectrum with the correct scan number */
    Spectrum* spectrum = collection->getSpectrum(scan_number);

    if( spectrum == NULL ){
      carp(CARP_WARNING, "Could not find scan number %i", scan_number);
      continue;
    }

    /* Print either the spectrum or stats. */
    if (!options){
      spectrum->print(stdout);

    } else {

      int charge_state_index = 0; 
      int charge_state_num = spectrum->getNumZStates();
      std::vector<SpectrumZState> zstates_array = spectrum->getZStates();
  
      printf("Scan number: %i\n", scan_number);
      printf("Precursor m/z:%.2f\n", spectrum->getPrecursorMz());
      printf("Total Ion Current:%.2f\n", spectrum->getTotalEnergy());
      printf("Base Peak Intensity:%.1f\n", 
             spectrum->getMaxPeakIntensity()); // base is max
      printf("Number of peaks:%d\n", spectrum->getNumPeaks());
      printf("Minimum m/z:%.1f\n", spectrum->getMinPeakMz());
      printf("Maximum m/z:%.1f\n", spectrum->getMaxPeakMz());
    
      for(charge_state_index=0; charge_state_index < charge_state_num; 
          ++charge_state_index){

        SpectrumZState& zstate = zstates_array[charge_state_index];
        FLOAT_T charged_mass = spectrum->getPrecursorMz() * (FLOAT_T)zstate.getCharge();

        printf("Charge state:%d\n", zstate.getCharge());
        printf("Neutral mass:%.2f\n", zstate.getNeutralMass());
        printf("Charged mass:%.2f\n", charged_mass);
        printf("M+H+ mass:%.2f\n", zstate.getSinglyChargedMass());
      }
    }
    delete spectrum;
    num_found++;
  }
  delete collection;

  carp(CARP_INFO, "Found %d spectra.\n", num_found);
  carp(CARP_INFO, "crux-get-ms2-spectrum finished.");
  return(0);
}
