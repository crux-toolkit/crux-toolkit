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
#include "unistd.h"

int main(int argc, char** argv){
  SPECTRUM_COLLECTION_T * collection; ///<spectrum collection
  SPECTRUM_T * spectrum; ///<the spectrum of interest
  BOOLEAN_T spectrum_found; ///<did we find the spectrum?
  int scan_num; ///< the query scan number
  FILE* output_file; ///< output file name
  BOOLEAN_T options = FALSE; ///< do we want options?
  char* USAGE = 
    "USAGE: get_ms2_spectrum [options]\n"
    "\t<scan number>\n"
    "\t<input filename>\n"
    "\t<output filename>\n";

  //check command line argument count
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
  
  //check if first argument is type int for scan number
  if (sscanf(argv[1], "%d", &scan_num) != 1) {
    fprintf (stderr, "ERROR:second arguement must be type int\n");
    fprintf (stderr, "%s", USAGE);
    exit(1);
  }

  //check if the input file exist
  if(access(argv[2], F_OK)){
    fprintf(stderr,"input file:\"%s\" could not be opened\n", argv[2]);
    exit(1);
  }
  
  //check if the output file doesn't aleady exist
  if(!access(argv[3], F_OK)){
    fprintf(stderr,"output file:\"%s\" already exist\n", argv[3]);
    exit(1);
  }

  //read input file
  collection = new_spectrum_collection(argv[2]);
  spectrum = allocate_spectrum();
  
  //search for spectrum with correct scan number
  if(get_spectrum_collection_spectrum(collection, scan_num, spectrum)){ 
    printf("Found a spectrum with correct scan number\n");
    output_file = fopen(argv[3], "w");
    print_spectrum(spectrum, output_file);
    fclose(output_file);
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

  int charge_state_index = 0; 
  int charge_state_num = get_spectrum_num_possible_z(spectrum);
  int* possible_z_array = get_spectrum_possible_z(spectrum);
  int possible_z;
  
  //print only if found spectrum and with option flag
  if(spectrum_found && options){
    printf("Precursor m/z:%.2f\n", get_spectrum_precursor_mz(spectrum));
    printf("Total Ion Current:%.2f\n", get_spectrum_total_energy(spectrum));
    printf("Base Peak Intensity:%.1f\n", get_spectrum_max_peak_intensity(spectrum)); //base is max
    printf("Number of peaks:%d\n", get_spectrum_num_peaks(spectrum));
    printf("Minimum m/z:%.1f\n", get_spectrum_min_peak_mz(spectrum));
    printf("Maximum m/z:%.1f\n", get_spectrum_max_peak_mz(spectrum));
    
    for(; charge_state_index < charge_state_num; ++charge_state_index){
      possible_z = possible_z_array[charge_state_index];
      printf("Charge state:%d\n", possible_z);
      printf("Neutral mass:%.2f\n", get_spectrum_neutral_mass(spectrum, possible_z));
      printf("Charged mass:%.2f\n", get_spectrum_mass(spectrum, possible_z));
      printf("M+H+ mass:%.2f\n", get_spectrum_singly_charged_mass(spectrum, possible_z));
    }
  }

  free(possible_z_array);
  free_spectrum(spectrum);
  free_spectrum_collection(collection);
  
  // failed to find spectrum
  if(!spectrum_found){
    return(1);
  }
  // exit success
  return(0);
}
