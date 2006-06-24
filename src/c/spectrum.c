/*****************************************************************************
 * \file spectrum.c
 * AUTHOR: Tobias P. Mann
 * CREATE DATE: 19 Sept 2003
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.4 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "spectrum.h"
#include "peak.h"
#include "utils.h"

/**
 * \define constants
 */
#define MAX_PEAKS 2000
#define MAX_CHARGE 2


/**
 * \struct spectrum 
 */
struct spectrum {
  int               first_scan;    ///< The number of the first scan
  int               last_scan;     ///< The number of the last scan
  int               id;            ///< A unique identifier
  SPECTRUM_TYPE_T   spectrum_type; ///< The type of spectrum.
  float             precursor_mz;  ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z;    ///< The possible charge states of this spectrum
  int               num_possible_z;///< The num of possible chare states of this spectrum
  PEAK_T*           peaks[MAX_PEAKS];         ///< The spectrum peaks
  float             min_peak_mz;   ///< The minimum m/z of all peaks
  float             max_peak_mz;   ///< The maximum m/z of all peaks
  int               num_peaks;     ///< The number of peaks
  double            total_energy;  ///< The sum of intensities in all peaks
  char*             filename;      ///< Optional filename
};    



BOOLEAN_T parse_S_line(SPECTRUM_T* spectrum, char* line, int buf_length);

BOOLEAN_T parse_Z_line(SPECTRUM_T* spectrum, char* line);

BOOLEAN_T add_possible_z(SPECTRUM_T* spectrum, int charge );

BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum, 
  float intensity,
  float location_mz );



/**
 * \returns An (empty) spectrum object.
 */
SPECTRUM_T* allocate_spectrum(){
  SPECTRUM_T* fresh_spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));
  fresh_spectrum->possible_z = (int*)mymalloc(sizeof(int) * MAX_CHARGE);
  return fresh_spectrum;
}

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */

/*
SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan
  int               last_scan,          ///< The number of the last scan
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum.
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z,         ///< The possible charge states of this spectrum
  char*             filename)          ///< Optional filename
{
  SPECTRUM_T* fresh_spectrum = allocate_spectrum();
  fresh_spectrum->first_scan = first_scan;
  fresh_spectrum->last_scan = last_scan;
  fresh_spectrum->spectrum_type = spectrum_type;
  fresh_spectrum->precursor_mz = precursor_mz;
  fresh_spectrum->possible_z = possible_z;
  fresh_spectrum->filename = filename; //optional, needs some work done on this.is it the file it reads?
} /////needs more work//////
*/

/**
 * Frees an allocated spectrum object.
 */
void free_spectrum (SPECTRUM_T* spectrum){
  int i = 0;
  free(spectrum->possible_z);
  for(; i < spectrum->num_peaks; ++i){
    free_peak(spectrum->peaks[i]); ///////think again!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  free(spectrum->filename);
  free(spectrum);
}


/**
 * Prints a spectrum object to file.
 */
void print_spectrum(SPECTRUM_T* spectrum, FILE* file){

}

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(SPECTRUM_T* spectrum){


}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T
 */
void copy_spectrum(
  SPECTRUM_T* src,
  SPECTRUM_T* dest)
{
  //copy each varible
  //copy each peak

}


/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 * 
 */
BOOLEAN_T parse_spectrum_file(
  SPECTRUM_T* spectrum,
  FILE* file)
{
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  float location_mz;
  float intensity;
  BOOLEAN_T record_S = FALSE; //check's if read S line
  BOOLEAN_T record_Z = FALSE;
  BOOLEAN_T start_add_peaks = FALSE;
  BOOLEAN_T file_format = FALSE;

  while( (line_length =  getline(&new_line, &buf_length, file)) != -1){
    if((!record_S || (record_S && start_add_peaks)) && 
       (new_line[0] == 'Z' ||  new_line[0] == 'I' ||
        new_line[0] == 'D' )){ // might add more checks
      file_format = FALSE;
      fprintf(stderr, "line: %s\n", new_line);
      break; //File format incorrect
    }
    else if(new_line[0] == 'S' && !record_S){
      record_S = TRUE;
      parse_S_line(spectrum, new_line, buf_length);
    }
    else if(new_line[0] == 'Z'){
      record_Z = TRUE;
      parse_Z_line(spectrum, new_line);
    }
    else if(new_line[0] == 'S' && start_add_peaks){ //start of next spectrum
      fprintf(stderr, "line: %s\n", new_line);
      break;
    }
    else if(record_Z && record_S &&  
            (sscanf(new_line,"%f %f", &location_mz, &intensity) == 2)){
      file_format = TRUE;
      start_add_peaks = TRUE;
      add_peak_to_spectrum(spectrum, intensity, location_mz);
    }
  }

  myfree(new_line);
  
  if(!file_format){ //File format incorrect
    fprintf(stderr, "incorrect file format\n");
    return FALSE;
  }
  return TRUE;
}

/**
 * Parses the 'S' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 * 
 */
BOOLEAN_T parse_S_line(SPECTRUM_T* spectrum, char* line, int buf_length){
  char spliced_line[buf_length];
  int line_index = 0;
  int spliced_line_index = 0;
  int first_scan;
  int last_scan;
  float precursor_mz;

  while(line[line_index] == 'S' || line[line_index] == '\t' ||
        line[line_index] == ' ' || line[line_index] == '0'){
    ++line_index;
  }
  while(line[line_index] != ' ' && line[line_index] != '\t'){
    spliced_line[spliced_line_index] =  line[line_index];
    ++spliced_line_index;
    ++line_index;
  }
  spliced_line[spliced_line_index] =  line[line_index];
  ++spliced_line_index;
  ++line_index;
  while(line[line_index] == '\t' || line[line_index] == ' ' || 
        line[line_index] == '0'){
    ++line_index;
  }
  while(line[line_index] !='\0'){
    spliced_line[spliced_line_index] =  line[line_index];
    ++spliced_line_index;
    ++line_index;
  }
  spliced_line[spliced_line_index] = '\0';
  if ( sscanf(spliced_line,"%i %i %f", 
              &first_scan, &last_scan, &precursor_mz) != 3) {
    fprintf(stderr,"Failed to parse S line %s",line);
    return FALSE;
  }
  //set_spectrum_first_scan( spectrum, first_scan);
  //set_spectrum_last_scan( spectrum, last_scan);
  //set_spectrum_precursor_mz( spectrum, precursor_mz);
  
  return TRUE;
}

/**
 * Parses the 'Z' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 * 
 */
BOOLEAN_T parse_Z_line(SPECTRUM_T* spectrum, char* line){
  int tokens;
  char line_name;
  int charge;
  float m_h_plus;
  if( (tokens = 
       sscanf(line, "%c %d %f", &line_name, &charge, &m_h_plus)) != 3){
    return FALSE;
  }  
  return add_possible_z(spectrum, charge);
}

/**
 * Adds a possible charge(z) to the spectrum.
 */
BOOLEAN_T add_possible_z(SPECTRUM_T* spectrum, int charge){
   int* possible_charge = (int *)mymalloc(sizeof(int));
   *possible_charge = charge;
   if(spectrum->num_possible_z < MAX_CHARGE){ // change to dynamic sometime...
     spectrum->possible_z[spectrum->num_possible_z] = *possible_charge; //make sure num_possible is initialized
     ++spectrum->num_possible_z;
     return TRUE;
   }
   return FALSE;
}


BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum, 
  float intensity, 
  float location_mz ) //return value???
{
  PEAK_T* peak = new_peak(intensity, location_mz);
  return add_peak(spectrum, peak); //think about return value on this..
}

/**
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_spectrum(
  SPECTRUM_T* spectrum,
  char*      filename)
{
  FILE* file;
  if ((file = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"File %s could not be opened",filename);
    return (FALSE);  //exit(1);
  }
  // might check if spectrum is NULL?? Close file??
  if(parse_spectrum_file(spectrum, file)){
    fclose(file);
    return TRUE;
  }
  return FALSE;
}

/**
 * Adds a peak to the spectrum.
 * update num_peaks, min_peak_mz, max_peak_mz, total_energy fields in spectrum
 */
BOOLEAN_T add_peak(SPECTRUM_T* spectrum, PEAK_T* peak){
  float location = peak_location(peak);
  
  if(spectrum->num_peaks < MAX_PEAKS){  // need to change it so that  it's dynamic
    spectrum->peaks[spectrum->num_peaks] = peak;
    ++spectrum->num_peaks;
  
    if(spectrum->num_peaks == 1 || 
       spectrum->min_peak_mz > location){
      spectrum->min_peak_mz = location;
    }
    if(spectrum->num_peaks == 1 || 
       spectrum->max_peak_mz < location){
      spectrum->max_peak_mz = location;
    }
    spectrum->total_energy += peak_intensity(peak);
    return TRUE;
  }
  return FALSE;
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

