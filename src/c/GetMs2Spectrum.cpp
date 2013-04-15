/**
 * \file GetMs2Spectrum.cpp
 *
 * AUTHOR: Manijeh Naser
 * CREATE DATE: January 31, 2012
 *
 * DESCRIPTION:brief searches a given ms2 file for the spectrum with the given
 * scan number.
 */

#include "GetMs2Spectrum.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <vector>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "parameter.h"
#include "carp.h"
#include "Spectrum.h"
#include "Peak.h"
#include "SpectrumCollectionFactory.h"
#include "WinCrux.h"
#include "ModifiedPeptidesIterator.h"

using namespace std;
using namespace Crux;

/**
 * \returns A blank GetMs2Spectrum object.
 */
GetMs2Spectrum::GetMs2Spectrum() {

}

/**
 * Destructorm
 */
GetMs2Spectrum::~GetMs2Spectrum() {
}

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
int GetMs2Spectrum :: main(int argc, char** argv){


  const char* option_list[] = { 
    "version", 
    "stats", 
    "verbosity",
    "spectrum-parser"};

 int num_options = sizeof(option_list) / sizeof(char*);

  
  const char* argument_list[] = {
    "scan number",
    "ms2 file"
  };
  int num_arguments = sizeof(argument_list) / sizeof(char*);
  initialize(argument_list,
	     num_arguments,
	     option_list,
	     num_options,
	     argc, argv);

  
  /* Get arguments */
  int min_scan;
  int max_scan;
  parse_scan_numbers(get_string_parameter("scan number"), &min_scan, &max_scan);
  fprintf(stderr, "Scanning from %d to %d.\n", min_scan, max_scan);

  const char* ms2_filename = get_string_parameter_pointer("ms2 file");
  carp(CARP_DETAILED_DEBUG, "ms2_filename: %s", ms2_filename);

  /* Get options */
  bool options = get_boolean_parameter("stats");

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
  
  return(0);
}
/**
 * \returns The command name for GetMs2SPectrum. 
 */
string GetMs2Spectrum::getName() {
  return "get-ms2-spectrum";
}

/**
 * \returns The description for GetMs2Spectrum.
 */
string GetMs2Spectrum::getDescription() {
  return "brief searches a given ms2 file for the spectrum with the given scan number.";
}

/**
 * \returns The enum of the application, default GET_MS2_SPECTRUM_COMMAND.
 */
COMMAND_T GetMs2Spectrum::getCommand() {
  return GET_MS2_SPECTRUM_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

