/*****************************************************************************
 * \file spectrum.c
 * AUTHOR: Chris Park
 * CREATE DATE:  June 22 2006
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.19 $
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
#include "mass.h"

/**
 * \define constants
 */
#define MAX_PEAKS 4000
#define MAX_CHARGE 2


/**
 * \struct spectrum 
 *  FIXME spectrum_type needs to be written..look at the file name .ms2
 */
struct spectrum {
  int               first_scan;    ///< The number of the first scan
  int               last_scan;     ///< The number of the last scan
  int               id;            ///< A unique identifier
  SPECTRUM_TYPE_T   spectrum_type; ///< The type of spectrum. 
  float             precursor_mz;  ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z;    ///< The possible charge states of this spectrum
  int               num_possible_z;///< The number of possible charge states of this spectrum
  PEAK_T*           peaks[MAX_PEAKS];         ///< The spectrum peaks
  float             min_peak_mz;   ///< The minimum m/z of all peaks
  float             max_peak_mz;   ///< The maximum m/z of all peaks
  int               num_peaks;     ///< The number of peaks
  double            total_energy;  ///< The sum of intensities in all peaks
  char*             filename;      ///< Optional filename
};    

/**
 * \struct peak_iterator
 */
struct peak_iterator {
  SPECTRUM_T* spectrum; ///< The spectrum whose peaks to iterate over. 
  int  peak_index;        ///< The index of the current peak
};

/**
 * Parses the 'S' line of a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_S_line(SPECTRUM_T* spectrum, char* line, int buf_length);

/**
 * Parses the 'Z' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_Z_line(SPECTRUM_T* spectrum, char* line);

/**
 * Adds a possible charge(z) to the spectrum.
 * Must not exceed the MAX_CHARGE capacity
 */
BOOLEAN_T add_possible_z(SPECTRUM_T* spectrum, int charge );

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

SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan
  int               last_scan,          ///< The number of the last scan
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum.
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z,         ///< The possible charge states of this spectrum
  int               num_possible_z,     ///< The number of possible charge states of this spectrum  
  char*             filename)          ///< Optional filename
{
  SPECTRUM_T* fresh_spectrum = allocate_spectrum();
  fresh_spectrum->first_scan = first_scan;
  fresh_spectrum->last_scan = last_scan;
  fresh_spectrum->spectrum_type = spectrum_type;
  fresh_spectrum->precursor_mz = precursor_mz;
  set_spectrum_new_possible_z(fresh_spectrum, possible_z, num_possible_z);
  set_spectrum_new_filename(fresh_spectrum, filename);
  return fresh_spectrum;
}


/**
 * Frees an allocated spectrum object.
 */
void free_spectrum (SPECTRUM_T* spectrum){
  int num_peaks_index = 0;
  free(spectrum->possible_z);
  for(; num_peaks_index < spectrum->num_peaks; ++num_peaks_index){
    free_peak(spectrum->peaks[num_peaks_index]);
  }
  free(spectrum->filename);
  free(spectrum);
}


/**
 * Prints a spectrum object to file.
 */
void print_spectrum(SPECTRUM_T* spectrum, FILE* file){
  int num_z_index = 0;
  int num_peak_index = 0;
  
  fprintf(file, "Filename: %s\n", spectrum->filename);
  fprintf(file, "S\t%06d\t%06d\t%.2f\n", 
         spectrum->first_scan,
         spectrum->last_scan,
         spectrum->precursor_mz);
  for(; num_z_index < spectrum->num_possible_z; ++num_z_index){
    fprintf(file, "Z\t%d\t%.2f\n", spectrum->possible_z[num_z_index],
            get_spectrum_singly_charged_mass(spectrum,
                                             spectrum->possible_z[num_z_index]));
  }
  for(; num_peak_index < spectrum->num_peaks; ++num_peak_index){
    fprintf(file, "%.1f %.1f\n", 
           peak_location(spectrum->peaks[num_peak_index]),
           peak_intensity(spectrum->peaks[num_peak_index]));
  }
}

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(SPECTRUM_T* spectrum){
  print_spectrum(spectrum, stdout);
}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T* dest
 */
void copy_spectrum(
  SPECTRUM_T* src,
  SPECTRUM_T* dest)
{
  int num_peak_index = 0;
  int* possible_z;
  char* new_filename;

  //copy each varible
  set_spectrum_first_scan(dest,get_spectrum_first_scan(src));
  set_spectrum_last_scan(dest,get_spectrum_last_scan(src));
  set_spectrum_id(dest,get_spectrum_id(src));
  set_spectrum_spectrum_type(dest,get_spectrum_spectrum_type(src));
  set_spectrum_precursor_mz(dest,get_spectrum_precursor_mz(src));
  
  //copy possible_z
  possible_z = get_spectrum_possible_z(src);
  set_spectrum_possible_z(dest,possible_z, 
                          get_spectrum_num_possible_z(src));
  free(possible_z);
  
  //copy filename
  new_filename = get_spectrum_filename(src);
  set_spectrum_filename(dest, new_filename);
  free(new_filename);
  
  //copy each peak
  for(; num_peak_index < get_spectrum_num_peaks(src); ++num_peak_index){
    add_peak_to_spectrum(dest, peak_intensity(src->peaks[num_peak_index]),
                         peak_location(src->peaks[num_peak_index])); 
  }
}


/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 * Must not include the Header line "H"
 * FIXME if need to read 'H' header line
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
  BOOLEAN_T record_S = FALSE; //check's if it read S line
  BOOLEAN_T record_Z = FALSE; //check's if it read Z line
  BOOLEAN_T start_add_peaks = FALSE; //check's if it started reading peaks
  BOOLEAN_T file_format = FALSE; //is the file format correct so far
  long file_index = ftell(file); //stores the location of the current working line in the file
  float test_float;
  char test_char;
  
  while( (line_length =  getline(&new_line, &buf_length, file)) != -1){
    // checks if 'S' is not the first line
    if((!record_S || (record_S && start_add_peaks)) && 
       (new_line[0] == 'Z' ||  
        new_line[0] == 'I' ||
        new_line[0] == 'D' )){
      file_format = FALSE;
      fprintf(stderr, "Incorrect order of line (S,Z, Peaks)\n");
      fprintf(stderr, "At line: %s", new_line);
      break; //File format incorrect
    }
    // Reads the 'S' line
    else if(new_line[0] == 'S' && !record_S){
      record_S = TRUE;
      if(!parse_S_line(spectrum, new_line, buf_length)){
        file_format = FALSE;
        break; //File format incorrect
      }
    }
    // Reads the 'Z' line 
    else if(new_line[0] == 'Z'){
      record_Z = TRUE;
      if(!parse_Z_line(spectrum, new_line)){
        file_format = FALSE;
        break; //File format incorrect
      }
    }
    // Stops, when encounters the start of next spectrum
    else if(new_line[0] == 'S' && start_add_peaks){ //start of next spectrum
      //fprintf(stderr, "line: %s\n", new_line);
      break;
    }

    //*****parse peak line******
    else if(new_line[0] != 'Z' &&  
            new_line[0] != 'I' &&
            new_line[0] != 'D' &&
            new_line[0] != '\n')
      {
        // checks if the peaks are in correct order of lines
        if((!record_Z || !record_S)){
          file_format = FALSE;
          fprintf(stderr, "Incorrect order of line (S,Z, Peaks)\n");
          fprintf(stderr, "At line: %s", new_line);
          break; //File format incorrect
        }
        // check for peak line format
        else if((sscanf(new_line,"%f %f %f",//test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_float) > 2)||
                (sscanf(new_line,"%f %f %c",//test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_char) > 2)||
                (sscanf(new_line,"%f %f",//test format:peak line has less than 2 fields
                        &test_float, &test_float) != 2)){
          file_format = FALSE;
          fprintf(stderr, "Incorrect peak line\n");
          fprintf(stderr, "At line: %s", new_line);
          break; //File format incorrect
        }
        // Reads the 'peak' lines, only if 'Z','S' line has been read
        else if(record_Z && record_S &&
                (sscanf(new_line,"%f %f", &location_mz, &intensity) == 2)){
          file_format = TRUE;
          start_add_peaks = TRUE;
          add_peak_to_spectrum(spectrum, intensity, location_mz);
        }
      }
    //*************************

    file_index = ftell(file); // updates the current working line location
  }
  
  // set the file pointer back to the start of the next 's' line
  fseek(file, file_index, SEEK_SET);
  myfree(new_line);
  set_spectrum_new_filename(spectrum,"test.ms2"); //FIXME temp field assignment

  //File format incorrect
  if(!file_format){ 
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
  float test_float;
  char test_char;

  //deletes empty space & 0
  while((line[line_index] !='\0') && 
        (line[line_index] == 'S' || 
         line[line_index] == '\t'||
         line[line_index] == ' ' || 
         line[line_index] == '0')){
    ++line_index;
  }
  // reads in line value
  while(line[line_index] !='\0' && 
        line[line_index] != ' ' && 
        line[line_index] != '\t'){
    spliced_line[spliced_line_index] =  line[line_index];
    ++spliced_line_index;
    ++line_index;
  }
  spliced_line[spliced_line_index] =  line[line_index];
  ++spliced_line_index;
  ++line_index;
  //deletes empty space & zeros
  while((line[line_index] !='\0') && 
        (line[line_index] == '\t' || 
         line[line_index] == ' ' || 
         line[line_index] == '0')){
    ++line_index;
  }
  // read last scan & precursor m/z
  while(line[line_index] !='\0'){
    spliced_line[spliced_line_index] =  line[line_index];
    ++spliced_line_index;
    ++line_index;
  }
  spliced_line[spliced_line_index] = '\0';
  
  // check if S line is in correct format
  if ( (sscanf(spliced_line,"%f %f %f %f",//test format:S line has more than 3 fields
               &test_float, &test_float, &test_float, &test_float) > 3) ||
       (sscanf(spliced_line,"%f %f %f %c",//test format:S line has more than 3 fields 
               &test_float, &test_float, &test_float, &test_char) > 3) ||
       (sscanf(spliced_line,"%i %i %f", // S line is parsed here
              &first_scan, &last_scan, &precursor_mz) != 3)) {
    fprintf(stderr,"Failed to parse 'S' line:\n %s",line);
    return FALSE;
  }
  set_spectrum_first_scan( spectrum, first_scan);
  set_spectrum_last_scan( spectrum, last_scan);
  set_spectrum_precursor_mz( spectrum, precursor_mz);
  
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
  float test_float;
  char test_char;

  // check if Z line is in correct format
  if( ((tokens =  // test format: Z line has less than 3 fields
        sscanf(line, "%c %f %f", &test_char, &test_float, &test_float)) < 3) ||
      ((tokens =   //test format: Z line has more than 3 fields
        sscanf(line, "%c %f %f %f", &test_char, &test_float, &test_float, &test_float)) >  3) ||
      ((tokens =  // test format: Z line has more than 3 fields
        sscanf(line, "%c %f %f %c", &test_char, &test_float, &test_float, &test_char)) >  3) ||
      (tokens = // Z line is parsed here
       sscanf(line, "%c %d %f", &line_name, &charge, &m_h_plus)) != 3){
    fprintf(stderr,"Failed to parse 'Z' line:\n %s",line);
    return FALSE;
  }  

  return add_possible_z(spectrum, charge);
}

/**
 * Adds a possible charge(z) to the spectrum.
 * Must not exceed the MAX_CHARGE capacity
 */
BOOLEAN_T add_possible_z(SPECTRUM_T* spectrum, int charge){
  int* possible_charge = (int *)mymalloc(sizeof(int));
   *possible_charge = charge;
   if(spectrum->num_possible_z < MAX_CHARGE){ // change to dynamic sometime...
     spectrum->possible_z[spectrum->num_possible_z] = *possible_charge; 
     ++spectrum->num_possible_z;
     free(possible_charge);
     return TRUE;
   }
   free(possible_charge);
   return FALSE;
}

/**
 * Adds a peak to the spectrum given a intensity and location
 * calls add_peak function
 */
BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum, 
  float intensity, 
  float location_mz )
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
  fclose(file);
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
    // is new peak the smallest peak
    if(spectrum->num_peaks == 1 || 
       spectrum->min_peak_mz > location){
      spectrum->min_peak_mz = location;
    }
    // is new peak the largest peak
    if(spectrum->num_peaks == 1 || 
       spectrum->max_peak_mz < location){
      spectrum->max_peak_mz = location;
    }
    //update total_energy
    spectrum->total_energy += peak_intensity(peak);
    return TRUE;
  }
  return FALSE;
}


/**
 * \returns the number of the first scan
 */
int get_spectrum_first_scan(SPECTRUM_T* spectrum){
  return spectrum->first_scan;
}

/**
 * \sets the number of the first scan
 */
void set_spectrum_first_scan(SPECTRUM_T* spectrum, int first_scan){
  spectrum->first_scan = first_scan;
}

/**
 * \returns the number of the last scan
 */
int get_spectrum_last_scan(SPECTRUM_T* spectrum){
  return spectrum->last_scan;
}

/**
 * \sets the number of the last scan
 */
void set_spectrum_last_scan(SPECTRUM_T* spectrum, int last_scan){
  spectrum->last_scan = last_scan;
}

/**
 * \returns the spectrum_id
 */
int get_spectrum_id(SPECTRUM_T* spectrum){
  return spectrum->id;
}

/**
 * \sets the spectrum_id
 */
void set_spectrum_id(SPECTRUM_T* spectrum, int id){
  spectrum->id = id;
}

/**
 * \returns the spectrum_type
 */
SPECTRUM_TYPE_T get_spectrum_spectrum_type(SPECTRUM_T* spectrum){
  return spectrum->spectrum_type;
}

/**
 * \sets the spectrum_type
 */
void set_spectrum_spectrum_type(SPECTRUM_T* spectrum, SPECTRUM_TYPE_T spectrum_type){
  spectrum->spectrum_type = spectrum_type;
}

/**
 * \returns the m/z of the precursor
 */
float get_spectrum_precursor_mz(SPECTRUM_T* spectrum){
  return spectrum->precursor_mz;
}

/**
 * \sets the m/z of the precursor
 */
void set_spectrum_precursor_mz(SPECTRUM_T* spectrum, float precursor_mz){
  spectrum->precursor_mz = precursor_mz;
}

/**
 * \returns the minimum m/z of all peaks
 */
float get_spectrum_min_peak_mz(SPECTRUM_T* spectrum){
  return spectrum->min_peak_mz;
}

/**
 * \returns the maximum m/z of all peaks
 */
float get_spectrum_max_peak_mz(SPECTRUM_T* spectrum){
  return spectrum->max_peak_mz;
}

/**
 * \returns the number of peaks
 */
int get_spectrum_num_peaks(SPECTRUM_T* spectrum){
  return spectrum->num_peaks;
}

/**
 * \returns the sum of intensities in all peaks
 */
double get_spectrum_total_energy(SPECTRUM_T* spectrum){
  return spectrum->total_energy;
}


/**
 * \returns the filename of the ms2 file the spectrum was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* get_spectrum_filename(SPECTRUM_T* spectrum){
  
  int filename_length = strlen(spectrum->filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename,spectrum->filename,filename_length);  
}

/**
 * \sets the filename of the ms2 file the spectrum was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void set_spectrum_filename(SPECTRUM_T* spectrum, char* filename){
  free(spectrum->filename);
  set_spectrum_new_filename(spectrum, filename);
}

/**  //////TESTME////
 * \sets the filename of the ms2 file the spectrum was parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_new_filename(SPECTRUM_T* spectrum, char* filename){
  int filename_length = strlen(filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  spectrum->filename =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * \returns the possible charge states of this spectrum
 * returns an int* to a heap allocated copy of the src spectrum
 * thus, the user must free the memory
 * number of possible charge states can be gained by 
 * the get_num_possible_z function.
 */
int* get_spectrum_possible_z(SPECTRUM_T* spectrum){
  int num_possible_z_index = 0;
  int* new_possible_z = 
    (int*)mymalloc(sizeof(int)*spectrum->num_possible_z);
  
  for(; num_possible_z_index < spectrum->num_possible_z; 
      ++num_possible_z_index){
  
    new_possible_z[num_possible_z_index]
      = spectrum->possible_z[num_possible_z_index];
  }
  return new_possible_z;
}
 
/**
 * \sets the possible charge states of this spectrum
 * this function should only be used when possible_z is set to NULL
 * to change existing possible_z use set_spectrum_possible_z()
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * updates the number of possible charge states field
 */
void set_spectrum_new_possible_z(SPECTRUM_T* spectrum, 
                                 int* possible_z, 
                                 int num_possible_z){
  int possible_z_index = 0;
  int* new_possible_z = 
    (int*)mymalloc(sizeof(int)*num_possible_z);
  
  for(; possible_z_index < num_possible_z; ++possible_z_index){
    new_possible_z[possible_z_index] = possible_z[possible_z_index];
  }
  
  spectrum->possible_z = new_possible_z;
  spectrum->num_possible_z = num_possible_z;

}

/**
 * \sets the possible charge states of this spectrum
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * frees the memory of the possible_z that is replaced
 * updates the number of possible charge states field
 */
void set_spectrum_possible_z(SPECTRUM_T* spectrum, 
                             int* possible_z, 
                             int num_possible_z){
  free(spectrum->possible_z);
  set_spectrum_new_possible_z(spectrum, possible_z, num_possible_z);
}

/**
 * \returns the number of possible charge states of this spectrum
 */
int get_spectrum_num_possible_z(SPECTRUM_T* spectrum){
  return spectrum->num_possible_z;
}

/**
 * \returns The intensity of the peak with the maximum intensity.
 */
float get_spectrum_max_peak_intensity(SPECTRUM_T* spectrum){
  int num_peak_index = 0;
  float max_intensity = -1;

  for(; num_peak_index < get_spectrum_num_peaks(spectrum); ++num_peak_index){
    if(max_intensity <= peak_intensity(spectrum->peaks[num_peak_index])){
      max_intensity = peak_intensity(spectrum->peaks[num_peak_index]);
    }
  }
  return max_intensity; 
}

/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
float get_spectrum_mass(SPECTRUM_T* spectrum, int charge){
  return get_spectrum_precursor_mz(spectrum)*charge;
}

/**
 * \returns The mass of the neutral precursor ion, according to the formula 
 * mass = m/z * charge - mass_H * charge
 */
float get_spectrum_neutral_mass(SPECTRUM_T* spectrum, int charge){
  return (get_spectrum_mass(spectrum, charge) - MASS_H*charge); //TESTME
}

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
float get_spectrum_singly_charged_mass(SPECTRUM_T* spectrum, int charge){
  return (get_spectrum_mass(spectrum, charge) - MASS_H*(charge-1));  //TESTME
}


/******************************************************************************/
// Iterator 

/**
 * Instantiates a new peak_iterator from a spectrum.
 * \returns a PEAK_ITERATOR_T object.
 */
PEAK_ITERATOR_T* new_peak_iterator(SPECTRUM_T* spectrum){
  PEAK_ITERATOR_T* peak_iterator = (PEAK_ITERATOR_T*)mycalloc(1,sizeof(PEAK_ITERATOR_T));
  peak_iterator->spectrum = spectrum;
  return peak_iterator;
}        

/**
 * Frees an allocated peak_iterator object.
 */
void free_peak_iterator(PEAK_ITERATOR_T* peak_iterator){
  free(peak_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peaks to iterate over, FALSE if not.
 */
BOOLEAN_T peak_iterator_has_next(PEAK_ITERATOR_T* peak_iterator){
  return  (peak_iterator->peak_index < get_spectrum_num_peaks(peak_iterator->spectrum));
}

/**
 * \returns The next peak object in the spectrum.
 */
PEAK_T* peak_iterator_next(PEAK_ITERATOR_T* peak_iterator){
  PEAK_T* next_peak = peak_iterator->spectrum->peaks[peak_iterator->peak_index];
  ++peak_iterator->peak_index;
  return next_peak;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

