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
#include "spectrum.h"
#include "peak.h"
#include "spectrum_collection.h"
#include "unistd.h"

int main(int argc, char** argv){
  SPECTRUM_COLLECTION_T * collection; ///<spectrum collection
  SPECTRUM_T * spectrum; ///<the spectrum of interest
  char* USAGE = 
    "USAGE: transform_ms2_spectra"
    "\t<input filename>\n"
    "SYNOPSIS: transform_ms2_spectra input_file > output_file\n";
  
  // check command line argument count
  if (argc != 2) {
    carp(CARP_FATAL, "%s", USAGE);
  }
  
  // check if the input file exist
  if(access(argv[1], F_OK)){
    carp(CARP_FATAL,"input file:\"%s\" could not be opened\n", argv[1]);
  }
  
  // read input file
  collection = new_spectrum_collection(argv[1]);
  if(!parse_spectrum_collection(collection)){
    carp(CARP_FATAL, "Failed to parse ms2 file: %s", ms2_file);
  }
  SPECTRUM_ITERATOR_T* iterator = new_spectrum_iterator(collection);
  

  int length = 1000000;
  int* seen = (int*) malloc(sizeof(int) * length);
  int idx;
  for (idx =0 ; idx < length; idx++){
    seen[idx] = 0;
  }

  while(spectrum_iterator_has_next(iterator)){
    spectrum = spectrum_iterator_next(iterator);
    int scan = get_spectrum_first_scan(spectrum);
    int* zs = get_spectrum_possible_z_pointer(spectrum);
    if (zs[0] == 1){
      add_possible_z(spectrum, 2);
      add_possible_z(spectrum, 3);
    } else if (zs[0] == 2){
      add_possible_z(spectrum, 1);
      add_possible_z(spectrum, 3);
    } else if (zs[0] == 3){
      add_possible_z(spectrum, 1);
      add_possible_z(spectrum, 2);
    }
    if (seen[scan] == 0){
      print_spectrum(spectrum, stdout);
      seen[scan] = 1;
    }
  }

  return(0);
}
