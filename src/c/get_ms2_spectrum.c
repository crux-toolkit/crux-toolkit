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
#include "spectrum.h"
#include "peak.h"
#include "spectrum_collection.h"

int main(int argc, char** argv){
  SPECTRUM_COLLECTION_T * collection; ///<spectrum collection
  SPECTRUM_T * spectrum; ///<the spectrum of interest
  BOOLEAN_T spectrum_found; ///<did we find the spectrum?
  int scan_num; ///< the query scan number
  char* USAGE = 
    "USAGE: get_ms2_spectrum\n"
    "\t<scan number>\n"
    "\t<input filename>\n"
    "\t<output filename>\n";

  // Parse command line options.
  if (argc != 4) {
    fprintf (stderr, "%s", USAGE);
    exit(1);
  }
  //FIXME
  // test argument type correct
  // file exist and doesn't exit..
  if (sscanf(argv[1], "%d", &scan_num) != 1) {
    fprintf (stderr, "ERROR:second arguement must be type int\n");
    fprintf (stderr, "%s", USAGE);
    exit(1);
  }

  //read input file
  collection = new_spectrum_collection(argv[2]);
  spectrum = allocate_spectrum();
  
  //search for spectrum with correct scan number
  if(get_spectrum_collection_spectrum(collection, scan_num, spectrum)){ 
    printf("Found a spectrum with correct scan number\n");
    print_spectrum_stdout(spectrum); //FIXME this is just for test
    spectrum_found = TRUE;
  }
  else{
    fprintf(stderr,
            "Warning:\n" 
            "The ms2 file %s does not contain the spectrum with scan number %d.\n", 
            argv[2], 
            scan_num);
    spectrum_found = FALSE;
  }
  

  free_spectrum(spectrum);
  free_spectrum_collection(collection);
  
  // failed to find spectrum
  if(!spectrum_found){
    return(1);
  }
  // exit success
  return(0);
}
