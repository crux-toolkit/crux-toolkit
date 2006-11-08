/*****************************************************************************
 * \file spectrum.c
 * AUTHOR: Chris Park
 * CREATE DATE:  June 22 2006
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.32 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "objects.h"
#include "spectrum.h"
#include "peak.h"
#include "utils.h"
#include "mass.h"
#include "parameter.h"
#include "scorer.h"

/**
 * \define constants
 */
#define MAX_PEAKS 4000
#define MAX_CHARGE 2

/**
 * \struct spectrum 
 * \brief A mass spectrum

 * A mass spectrum consists mainly of a list of peak objects along with
 * some identifying information. A single spectrum is generated from one 
 * or more "scans" of the mass spectrometer; each scan is identified by 
 * a unique increasing positive integer. The range of scans that
 * generated a particular spectrum are indicated by the member variables 
 * "first_scan" and "last_scan". In addition to scan information, 
 * a tandem fragmentation mass spectrum has information 
 * about the m/z of the intact ion that generated the spectrum, which is
 * indicated by "precursor_mz" member variable.
 * Also, while the m/z of particular spectrum is known, the charge state of
 * the originating ion is unknown; the possible charge states of the 
 * precursor ion is stored "possible_z" and "num_possible_z". 
 * Finally, some summary information that can be derived from the spectrum
 * peaks but is convenient to have is stored as "min_peak_mz",
 * "max_peak_mz", and "total_energy".
 */
struct spectrum {
  int               first_scan;    ///< The number of the first scan
  int               last_scan;     ///< The number of the last scan
  int               id;            ///< A unique identifier FIXME, this field is not set when parsing..
  SPECTRUM_TYPE_T   spectrum_type; ///< The type of spectrum. 
  float             precursor_mz;  ///< The m/z of the precursor (for MS-MS spectra)
  int*              possible_z;    ///< The possible charge states of this spectrum
  int               num_possible_z;///< The number of possible charge states of this spectrum
  PEAK_T*           peaks;         ///< The spectrum peaks
  float             min_peak_mz;   ///< The minimum m/z of all peaks
  float             max_peak_mz;   ///< The maximum m/z of all peaks
  int               num_peaks;     ///< The number of peaks
  double            total_energy;  ///< The sum of intensities in all peaks
  char*             filename;      ///< Optional filename

  //spectrum scoring fields, all these fields are valid only after make sum array
  //BOOLEAN_T         sum_array_exist; ///< does the sum array exist?
  //float*            intensity_sum_array; ///< intesity sum array used for 
  //float             sp_sum_resolution;  ///< the resolution of the intesity sum, the closure of which the array was created  
  float             max_intensity;   ///< the maximum intensity of peaks, only valid after sum array is created
  BOOLEAN_T         ready_for_sp; ///< has this spectrum been processed for SP?
};    

/**
 * \struct peak_iterator
 * \brief Object to iterate over the peaks in a spectrum.
 */
struct peak_iterator {
  SPECTRUM_T* spectrum; ///< The spectrum whose peaks to iterate over. 
  int  peak_index;        ///< The index of the current peak
};

/*
 * Parses the 'S' line of a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_S_line(
  SPECTRUM_T* spectrum, ///< place S line into this spectrum -out 
  char* line, ///< 'S' line to parse -in
  int buf_length ///< line length -in
);

/*
 * Parses the 'Z' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_Z_line(
  SPECTRUM_T* spectrum, ///< place Z line into this spectrum -out  
  char* line  ///< 'Z' line to parse -in
);

/*
 * Adds a possible charge(z) to the spectrum.
 * Must not exceed the MAX_CHARGE capacity
 */
BOOLEAN_T add_possible_z(
  SPECTRUM_T* spectrum,  ///< place Z line into this spectrum -out   
  int charge  ///< charge to add
);

/**
 * \returns An (empty) spectrum object.
 */
SPECTRUM_T* allocate_spectrum(void){
  SPECTRUM_T* fresh_spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));
  fresh_spectrum->possible_z = (int*)mymalloc(sizeof(int) * MAX_CHARGE);
  fresh_spectrum->peaks = allocate_peak_array(MAX_PEAKS);
  //fresh_spectrum->sum_array_exist = FALSE;
  fresh_spectrum->ready_for_sp = FALSE;
  return fresh_spectrum;
}

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */ 

SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan -in
  int               last_scan,          ///< The number of the last scan -in
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum. -in
  float             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra) -in
  int*              possible_z,         ///< The possible charge states of this spectrum  -in
  int               num_possible_z,     ///< The number of possible charge states of this spectrum  -in
  char*             filename)          ///< Optional filename -in
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
void free_spectrum (
  SPECTRUM_T* spectrum ///< the spectrum to free -in
)
{
  free(spectrum->possible_z);
  free(spectrum->filename);
  free(spectrum->peaks);
  //free(spectrum->intensity_sum_array);
  free(spectrum);
}


/**
 * Prints a spectrum object to file.
 */
void print_spectrum(
  SPECTRUM_T* spectrum, ///< spectrum to print -in
  FILE* file ///< output file to print at -out
  )
{
  int num_z_index = 0;
  int num_peak_index = 0;

  fprintf(file, "Filename: %s\n", spectrum->filename);
  fprintf(file, "S\t%06d\t%06d\t%.2f\n", 
         spectrum->first_scan,
         spectrum->last_scan,
         spectrum->precursor_mz);
  //print 'Z' line
  for(; num_z_index < spectrum->num_possible_z; ++num_z_index){
    fprintf(file, "Z\t%d\t%.2f\n", spectrum->possible_z[num_z_index],
            get_spectrum_singly_charged_mass(spectrum,
                                             spectrum->possible_z[num_z_index]));
  }
  //print peaks
  for(; num_peak_index < spectrum->num_peaks; ++num_peak_index){
    fprintf(file, "%.1f %.1f\n", 
            get_peak_location(find_peak(spectrum->peaks, num_peak_index)),
            get_peak_intensity(find_peak(spectrum->peaks, num_peak_index)));
  }
}

/**
 * Prints a spectrum object to STDOUT.
 */
void print_spectrum_stdout(
  SPECTRUM_T* spectrum ///< the spectrum to print -in
  )
{
  print_spectrum(spectrum, stdout);
}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T* dest
 * doesn't copy the sum array related fields
 */
void copy_spectrum(
  SPECTRUM_T* src, ///< the source spectrum -in
  SPECTRUM_T* dest ///< the destination spectrum -out
  )
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
    add_peak_to_spectrum(dest, get_peak_intensity(find_peak(src->peaks, num_peak_index)),
                         get_peak_location(find_peak(src->peaks, num_peak_index))); 
  }
}

/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 * Skips Header line "H"
 * FIXME if need to read 'H','I' header line, does not parse ID
 */
BOOLEAN_T parse_spectrum_file(
  SPECTRUM_T* spectrum, ///< spectrum to parse the information into -out
  FILE* file ///< the input file stream -in
  )
{
  long file_index = ftell(file); //stores the location of the current working line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  float location_mz;
  float intensity;
  BOOLEAN_T record_S = FALSE; //check's if it read S line
  BOOLEAN_T record_Z = FALSE; //check's if it read Z line
  BOOLEAN_T start_add_peaks = FALSE; //check's if it started reading peaks
  BOOLEAN_T file_format = FALSE; //is the file format correct so far
  
  float test_float;
  char test_char;
  
  while( (line_length =  getline(&new_line, &buf_length, file)) != -1){
    //skip header line
    //if(new_line[0] == 'H'){
    //  file_index = ftell(file);
    //  continue;
    //}
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
    // Stops, when encounters the start of next spectrum 'S' line
    else if(new_line[0] == 'S' && start_add_peaks){ //start of next spectrum
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
  
  //No more spectrum in .ms file
  if(!record_S && !file_format){
    return FALSE;
  }
  
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
BOOLEAN_T parse_S_line(
  SPECTRUM_T* spectrum, ///< place S line into this spectrum -out 
  char* line, ///< 'S' line to parse -in
  int buf_length ///< line length -in
  )
{
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
BOOLEAN_T parse_Z_line(
  SPECTRUM_T* spectrum, ///< place Z line into this spectrum -out  
  char* line  ///< 'Z' line to parse -in
  )
{
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
BOOLEAN_T add_possible_z(
  SPECTRUM_T* spectrum,  ///< place Z line into this spectrum -out   
  int charge  ///< charge to add
  )
{
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
 * calls update_spectrum_fields to update num_peaks, min_peak ...
 */
BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum,///< spectrum to add the peak to -out 
  float intensity, ///< the intensity of peak to add -in
  float location_mz ///< the location of peak to add -in
  )
{
  if(spectrum->num_peaks < MAX_PEAKS){  //FIXME someday change it to be dynamic
    set_peak_intensity(find_peak(spectrum->peaks, spectrum->num_peaks), intensity);
    set_peak_location(find_peak(spectrum->peaks, spectrum->num_peaks), location_mz);
    update_spectrum_fields(spectrum, intensity, location_mz);
    return TRUE;
  }

  return FALSE;
}

/**
 * \returns TRUE if success. FALSE if failure.
 */
BOOLEAN_T parse_spectrum(
  SPECTRUM_T* spectrum, ///< spectrum to parse the information into -out
  char*      filename ///< the file to parse -in
  )
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
 * updates num_peaks, min_peak_mz, max_peak_mz, total_energy fields in spectrum
 */
void update_spectrum_fields(
  SPECTRUM_T* spectrum, ///< the spectrum fields to update -out
  float intensity, ///< the intensity of the peak that has been added -in
  float location ///< the location of the peak that has been added -in
  )
{
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
  spectrum->total_energy += intensity;
}


/**
 * \returns the number of the first scan
 */
int get_spectrum_first_scan(
  SPECTRUM_T* spectrum ///< the spectrum to query the first scan -in
  )
{
  return spectrum->first_scan;
}

/**
 * sets the number of the first scan
 */
void set_spectrum_first_scan(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the first scan -out
  int first_scan ///< the first_scan -in
  )
{
  spectrum->first_scan = first_scan;
}

/**
 * \returns the number of the last scan
 */
int get_spectrum_last_scan(
  SPECTRUM_T* spectrum  ///< the spectrum to query the last scan -in
  )
{
  return spectrum->last_scan;
}

/**
 * sets the number of the last scan
 */
void set_spectrum_last_scan(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the last scan -out
  int last_scan ///< the last scan -in
  )
{
  spectrum->last_scan = last_scan;
}

/**
 * \returns the spectrum_id
 */
int get_spectrum_id(
  SPECTRUM_T* spectrum  ///< the spectrum to query the id -in
  )
{
  return spectrum->id;
}

/**
 * sets the spectrum_id
 */
void set_spectrum_id(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the id -out
  int id ///< the id -in
  )
{
  spectrum->id = id;
}

/**
 * \returns the spectrum_type
 */
SPECTRUM_TYPE_T get_spectrum_spectrum_type(
  SPECTRUM_T* spectrum  ///< the spectrum to query the spectrum_type -in
  )
{
  return spectrum->spectrum_type;
}

/**
 * sets the spectrum_type
 */
void set_spectrum_spectrum_type(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the spectrum_type -out
  SPECTRUM_TYPE_T spectrum_type ///< the spectrum type -in
  )
{
  spectrum->spectrum_type = spectrum_type;
}

/**
 * \returns the m/z of the precursor
 */
float get_spectrum_precursor_mz(
  SPECTRUM_T* spectrum  ///< the spectrum to query the precursor_mz -in
  )
{
  return spectrum->precursor_mz;
}

/**
 * sets the m/z of the precursor
 */
void set_spectrum_precursor_mz(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the precursor_mz -out
  float precursor_mz ///< the precursor_mz -in
  )
{
  spectrum->precursor_mz = precursor_mz;
}

/**
 * \returns the minimum m/z of all peaks
 */
float get_spectrum_min_peak_mz(
  SPECTRUM_T* spectrum ///< the spectrum to query min_peak_mz -in
  )
{
  return spectrum->min_peak_mz;
}

/**
 * \returns the maximum m/z of all peaks
 */
float get_spectrum_max_peak_mz(
  SPECTRUM_T* spectrum  ///< the spectrum to query max_peak_mz -in
  )
{
  return spectrum->max_peak_mz;
}

/**
 * \returns the number of peaks
 */
int get_spectrum_num_peaks(
  SPECTRUM_T* spectrum  ///< the spectrum to query number of peaks -in
  )
{
  return spectrum->num_peaks;
}

/**
 * \returns the sum of intensities in all peaks
 */
double get_spectrum_total_energy(
  SPECTRUM_T* spectrum  ///< the spectrum to query total energy -in
  )
{
  return spectrum->total_energy;
}


/**
 * \returns the filename of the ms2 file the spectrum was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* get_spectrum_filename(
  SPECTRUM_T* spectrum  ///< the spectrum to query filename -in
  )
{
  
  int filename_length = strlen(spectrum->filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename,spectrum->filename,filename_length);  
}

/**
 * sets the filename of the ms2 file the spectrum was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void set_spectrum_filename(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the filename -out
  char* filename ///< the filename -in
  )
{
  free(spectrum->filename);
  set_spectrum_new_filename(spectrum, filename);
}

/**
 * sets the filename of the ms2 file the spectrum was parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_new_filename(
  SPECTRUM_T* spectrum,   ///< the spectrum to set the filename -out
  char* filename ///< the filename -in
  )
{
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
int* get_spectrum_possible_z(
  SPECTRUM_T* spectrum  ///< the spectrum to query possible z -in
  )
{
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
 * sets the possible charge states of this spectrum
 * this function should only be used when possible_z is set to NULL
 * to change existing possible_z use set_spectrum_possible_z()
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * updates the number of possible charge states field
 */
void set_spectrum_new_possible_z(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the new_possible_z -out
  int* possible_z, ///< possible z array -in
  int num_possible_z ///< possible z array size -in
  )
{
  
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
 * sets the possible charge states of this spectrum
 * the function copies the possible_z into a heap allocated memory
 * num_possible_z must match the array size of possible_z 
 * frees the memory of the possible_z that is replaced
 * updates the number of possible charge states field
 */
void set_spectrum_possible_z(
  SPECTRUM_T* spectrum,  ///< the spectrum to set the new_possible_z -out
  int* possible_z, ///< possible z array -in
  int num_possible_z ///< possible z array size -in
  )
{
  free(spectrum->possible_z);
  set_spectrum_new_possible_z(spectrum, possible_z, num_possible_z);
}

/**
 * \returns the number of possible charge states of this spectrum
 */
int get_spectrum_num_possible_z(
  SPECTRUM_T* spectrum  ///< the spectrum to query length of possible_z array -in
  )
{
  return spectrum->num_possible_z;
}

/**
 * \returns The intensity of the peak with the maximum intensity.
 */
float get_spectrum_max_peak_intensity(
  SPECTRUM_T* spectrum  ///< the spectrum to query maximum peak intensity -in
  )
{
  int num_peak_index = 0;
  float max_intensity = -1;

  for(; num_peak_index < get_spectrum_num_peaks(spectrum); ++num_peak_index){
    if(max_intensity <= get_peak_intensity(find_peak(spectrum->peaks, num_peak_index))){
      max_intensity = get_peak_intensity(find_peak(spectrum->peaks, num_peak_index));
    }
  }
  return max_intensity; 
}


/**
 * Only should be used after process spectrum, other times use get_spectrum_max_peak_intensity
 * \returns The intensity of the peak with the maximum intensity.
 */
float get_spectrum_max_intensity(
  SPECTRUM_T* spectrum  ///< the spectrum to query maximum peak intensity -in
  )
{ 
  // check if max_intensity value is valid
  if(!spectrum->ready_for_sp){
    carp(CARP_ERROR, "cannot use max_intensity must first process spectrum");
    exit(1);
  }
  
  return spectrum->max_intensity;
}

/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
float get_spectrum_mass(
  SPECTRUM_T* spectrum,  ///< the spectrum to query spectrum mass -in
  int charge ///< the charge of precursor ion -in
  )
{
  return get_spectrum_precursor_mz(spectrum)*charge;
}

/**
 * \returns The mass of the neutral precursor ion, according to the formula 
 * mass = m/z * charge - mass_H * charge
 */
float get_spectrum_neutral_mass(
  SPECTRUM_T* spectrum,  ///< the spectrum to query neutral_mass -in
  int charge ///< the charge of precursor ion -in
  )
{
  return (get_spectrum_mass(spectrum, charge) - MASS_H*charge); //TESTME
}

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
float get_spectrum_singly_charged_mass(
  SPECTRUM_T* spectrum,  ///< the spectrum to query charged_mass -in
  int charge ///< the charge of the precursor ion -in
  )
{
  return (get_spectrum_mass(spectrum, charge) - MASS_H*(charge-1));  //TESTME
}



/******************************************************************************/
// Iterator 

/**
 * Instantiates a new peak_iterator from a spectrum.
 * \returns a PEAK_ITERATOR_T object.
 */
PEAK_ITERATOR_T* new_peak_iterator(
  SPECTRUM_T* spectrum ///< the spectrum peaks to iterate -in
  )
{
  PEAK_ITERATOR_T* peak_iterator = (PEAK_ITERATOR_T*)mycalloc(1,sizeof(PEAK_ITERATOR_T));
  peak_iterator->spectrum = spectrum;
  return peak_iterator;
}        

/**
 * Frees an allocated peak_iterator object.
 */
void free_peak_iterator(
  PEAK_ITERATOR_T* peak_iterator ///< the interator to free -in
  )
{
  free(peak_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peaks to iterate over, FALSE if not.
 */
BOOLEAN_T peak_iterator_has_next(
  PEAK_ITERATOR_T* peak_iterator ///< the interator for the peaks -in
  )
{
  return  (peak_iterator->peak_index < get_spectrum_num_peaks(peak_iterator->spectrum));
}

/**
 * \returns The next peak object in the spectrum.
 */
PEAK_T* peak_iterator_next(
  PEAK_ITERATOR_T* peak_iterator  ///< the interator for the peaks -in
  )
{
  PEAK_T* next_peak = find_peak(peak_iterator->spectrum->peaks, peak_iterator->peak_index);
  ++peak_iterator->peak_index;
  return next_peak;
}



/******************************************************************
 * The scoring relating methods
 ******************************************************************/


/**
 * For SP!
 * Adds a peak to the spectrum given a intensity and location, if
 * the peak has a unique m/z, if not, sets the peak with the same m/z with the
 * largest intensity. Thus, only one peak per m/z
 * In spectrum fields, only updates the num peaks, other fields are left unchanged
 */
BOOLEAN_T add_or_update_peak_to_spectrum_for_sp(
  SPECTRUM_T* spectrum,///< spectrum to add the peak to -out 
  float intensity, ///< the intensity of peak to add -in
  int location_mz ///< the location of peak to add -in
  )
{
  PEAK_T* prior_peak = NULL;
  
  if(spectrum->num_peaks < MAX_PEAKS){  //FIXME someday change it to be dynamic
    //check if only there are any prior peaks added
    if(spectrum->num_peaks > 0){
      prior_peak = find_peak(spectrum->peaks, spectrum->num_peaks-1);
      
      //there is a peak with same m/z
      if((int)get_peak_location(prior_peak) == location_mz){
        //if peak's intensity is smaller than the new intensity, replace with new intensity
        if(get_peak_intensity(prior_peak) < intensity){
          set_peak_intensity(prior_peak, intensity);
        }
        return TRUE;
      }
    }
        
    set_peak_intensity(find_peak(spectrum->peaks, spectrum->num_peaks), intensity);
    set_peak_location(find_peak(spectrum->peaks, spectrum->num_peaks), location_mz);
    ++spectrum->num_peaks;
    return TRUE;
  }
  
  return FALSE;
}

/**
 * Copies spectrum object src to dest.
 * must pass in a memory allocated SPECTRUM_T* dest
 * excludes peaks that are 15u around the precuros ion and larger than 50+experimental mass
 * Square roots all the intensity
 * sets the maximum_peak_intensity field, rounds all peak location to nearest int 
 */
void copy_spectrum_for_sp(
  SPECTRUM_T* src, ///< the source spectrum -in
  SPECTRUM_T* dest ///< the destination spectrum -out
  )
{
  int num_peak_index = 0;
  int* possible_z;
  char* new_filename;

  //FIXME...might not be correct!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  float experimental_mass_cut_off = get_spectrum_precursor_mz(src)*src->possible_z[0] + 50;

  //FIXME maybe add 0.5
  int precursor_mz = (int)(get_spectrum_precursor_mz(src)+0.5);
  float max_intensity = 0;
  float working_intensity = 0;
  float peak_location = 0;
  int working_location = 0;

  //copy necessary varibles
  set_spectrum_spectrum_type(dest,get_spectrum_spectrum_type(src));
  set_spectrum_precursor_mz(dest, precursor_mz);

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
    
    //get peak location
    peak_location = get_peak_location(find_peak(src->peaks, num_peak_index));

    //do not add any ion above the experimental mass
    if(peak_location > experimental_mass_cut_off){
      break;
    }
    
    //round the location to nearest int
    working_location = (int)(peak_location + 0.5);
    
    //skip all peaks within 15u of precursor_mz
    if(working_location <= precursor_mz + 15 && working_location >= precursor_mz - 15){
      continue;
    }
    
    //square root peak's intensity
    working_intensity = sqrt(get_peak_intensity(find_peak(src->peaks, num_peak_index)));

    //set maximum intensity if largest so far
    if(working_intensity > max_intensity){
      max_intensity = working_intensity;
    }
    
    //add or updates an existing peak to dest
    add_or_update_peak_to_spectrum_for_sp(dest, working_intensity, working_location); 
  }

  //set maximum peak intensity
  dest->max_intensity = max_intensity;
}


/**
 * process the spectrum, for SP
 */
void sp_process_spectrum(
  SPECTRUM_T* dest ///< the spectrum to processes -in/out
  )
{
  
  //keep only 200 most abundant ions
  if(dest->num_peaks > 200){
    dest->num_peaks = 200;
  }
  
  //re sort back to intensity order
  //sort the peak by location
  sort_peaks(dest->peaks, dest->num_peaks, _PEAK_LOCATION);    

}

/**
 * process the spectrum, according the score type
 *\returns a new spectrum that has been preprocessed
 */
SPECTRUM_T* process_spectrum(
  SPECTRUM_T* spectrum, ///< the spectrum to processes -in
  SCORER_TYPE_T score_type ///< the score type to which the spectrum should be sorted -in
  )
{
  //copy spectrum to dest, make a duplicate spectrum to produce a preprocessed spectrum
  SPECTRUM_T* dest = (SPECTRUM_T*)allocate_spectrum();
  
  if(score_type == SP){
    //copy the spectrum, excluding peaks 15u around the precursor ion
    copy_spectrum_for_sp(spectrum, dest);
    
    //sort the peak by intensity
    sort_peaks(dest->peaks, dest->num_peaks, _PEAK_INTENSITY);
    
    sp_process_spectrum(dest);
    
    //now ready for sp scoring
    dest->ready_for_sp = TRUE;
  }
  else if(score_type == XCORR){
    //add code here for xcore
  }
  
  //debug
  //print_spectrum_stdout(dest);

  return dest;
}

/******************************
 * OLD routines
 ******************************/

/**
 *\returns the index of which the intensity sum array corresponds to the theoretical mz
 */
/*
int hash_sum_array(
  SCORER_T* scorer,        ///< the scorer object -in  
  //SPECTRUM_T* spectrum, ///< the spectrum to query the intensity sum -in
  float peak_mz ///< the theoretical peaks mz -in
  )
{
  //set sum array bounderies, and other needed variables
  float min_mz = get_scorer_sp_min_mz(scorer);
  float max_mz = get_scorer_sp_max_mz(scorer);
  float sp_array_resolution = get_scorer_sp_array_resolution(scorer);
  
  //check if mz within mz handle limit
  if(peak_mz < min_mz || peak_mz > max_mz){
    carp(CARP_ERROR, "peaks mz not supported, bellow minimum mz");
    exit(1); //FIXME maybe should not die..
  }
  
  //get base index fraction
  float index_fraction = (peak_mz - min_mz)/sp_array_resolution;
  
  //get the base index
  int hash_index = (int)index_fraction;
  
  //if the intensity is closer to the next index, return next index
  
  //if(index_fraction - hash_index > sp_array_resolution){
  //  ++hash_index;
  // }
  
  return hash_index;
}
*/

/**
 * adds the intensity to the summ array
 * adds to all the bins that the interval of mz fits the intensity
 */
/*
void add_intensity_to_sum_array(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum, ///< the spectrum to add the intensity -in
  int array_size, ///< the size of the array -in
  float intensity, ///< the intensity to add -in
  float mz ///< the peaks mz whose intensity to add -in
  )
{
  float* sum_array = spectrum->intensity_sum_array; ///< the sum array to add 
  float sp_sum_resolution = spectrum->sp_sum_resolution; ///< the resolution of the intensity sum 
  float sp_array_resolution = get_scorer_sp_array_resolution(scorer);
  float min_mz = get_scorer_sp_min_mz(scorer);

  //get the base index
  int hash_index = hash_sum_array(scorer, mz);
  int sum_index = hash_index;

  sum_array[sum_index] += intensity;

  //add to bins after hashed bin
  while(mz <= (min_mz + (++sum_index+1)*sp_array_resolution) + sp_sum_resolution &&
        mz >= (min_mz + (sum_index)*sp_array_resolution) - sp_sum_resolution &&
        sum_index < array_size){
      sum_array[sum_index] += intensity;
  }

  sum_index = hash_index;       //check for index going over 0 or max

  //add to bins before hased bin
  while(mz <= (min_mz + (--sum_index + 1)*sp_array_resolution) + sp_sum_resolution &&
        mz >= (min_mz + (sum_index)*sp_array_resolution) - sp_sum_resolution &&
        sum_index >= 0){
    sum_array[sum_index] += intensity;
  }
}
*/


/**
 * process the spectrum, according the score type
 *\returns a new spectrum that has been preprocessed
 */
/*
SPECTRUM_T* process_spectrum(
  SPECTRUM_T* spectrum, ///< the spectrum to processes -in
  SCORER_TYPE_T score_type ///< the score type to which the spectrum should be sorted -in
  )
{
  //copy spectrum to dest, make a duplicate spectrum to produce a preprocessed spectrum
  SPECTRUM_T* dest = (SPECTRUM_T*)allocate_spectrum();

  //copy the spectrum, excluding peaks 10u around the precursor ion
  //also sets the maximum intensity field in dest
  copy_spectrum_for_process(spectrum, dest);
  
  //sort the peak by intensity
  sort_peaks(dest->peaks, dest->num_peaks, _PEAK_INTENSITY);
  
  if(score_type == SP){
    sp_process_spectrum(dest);
  }
  else if(score_type == XCORR){
    //add code here for xcore
  }
  
  //debug
  //print_spectrum_stdout(dest);

  return dest;
}
*/

/**
 * process the spectrum
 *\returns a new spectrum that has been preprocessed, frees old one.
 */
/*
SPECTRUM_T* process_spectrum(
  SPECTRUM_T* spectrum ///< the spectrum to processes -in
  )
{
  
  dest = (SPECTRUM_T*)allocate_spectrum();

  int peak_idx = 0;
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
  for(; num_peak_index < src->num_peaks; ++num_peak_index){
    peak = find_peak(spectrum->peaks, peak_idx);
    intensity = get_peak_intensity(peak);
    peak_mz = get_peak_location(peak);

    add_peak_to_spectrum(dest, get_peak_intensity(find_peak(src->peaks, num_peak_index)),
                         get_peak_location(find_peak(src->peaks, num_peak_index))); 
  }
  */
  /*
  int peak_idx = 0;
  float intensity = 0;
  PEAK_T* peak = NULL;
  PEAK_T* peak_around = NULL;
  //float max_intensity = 0;
  float peak_mz = 0;
  
  //iterate over all peaks and add the intensity to the correct array sum bins
  for(; peak_idx < spectrum->num_peaks; ++peak_idx){
    peak = find_peak(spectrum->peaks, peak_idx);
    intensity = get_peak_intensity(peak);
    peak_mz = get_peak_location(peak);

    peak_around = find_peak(spectrum->peaks, peak_idx+1);
    if(peak_idx != spectrum->num_peaks-1 && peak_mz+1 > get_peak_location(peak_around)){
      if(intensity > get_peak_intensity(peak_around)){
        set_peak_intensity(peak_around, intensity);
      }
      else{
        set_peak_intensity(peak, get_peak_intensity(peak_around));
      }
    }

    peak_around = find_peak(spectrum->peaks, peak_idx-1);
    if(peak_idx != 0 && peak_mz-1 < get_peak_location(peak_around)){
      if(intensity > get_peak_intensity(peak_around)){
        set_peak_intensity(peak_around, intensity);
      }
      else{
        set_peak_intensity(peak, get_peak_intensity(peak_around));
      }
    }
  }
}
*/

/**
 *
 *\set the peak intensity to the highest of it's +/- 1u peaks
 */
/*
void equalize_peak(
  SPECTRUM_T* spectrum, ///< the spectrum to query the intensity sum -in
  int peak_idx, ///< the current peak idx -in
  float* intensity, ///< the peak's intensity to equalize -in/out
  float peak_mz,  ///< the peak's location -in
  float sp_equalize_resolution ///< the eualization resolution -in
  )
{
  int working_idx = peak_idx+1;
  PEAK_T* peak = NULL;
  float mz = 0;
  float working_intensity;

  //FIXME, this is to avoid float rounding problems (ex, 0.9999)
  sp_equalize_resolution += 0.01;

  //iterate over all peaks and add the intensity to the correct array sum bins
  for(; working_idx < spectrum->num_peaks; ++working_idx){
    peak = find_peak(spectrum->peaks, working_idx);
    mz = get_peak_location(peak);
    
    //too small go to next peak
    if(mz < (peak_mz - sp_equalize_resolution)){
      continue;
    }
    
    //beyond range, break
    if(mz > (peak_mz + sp_equalize_resolution)){
      break;
    }

    working_intensity = get_peak_intensity(peak);
    
    //check it the intensity is larger than the peak to equalize
    if(*intensity < working_intensity){
      *intensity = working_intensity;
    }
  }

  working_idx = peak_idx-1;

  //iterate over all peaks and add the intensity to the correct array sum bins
  for(; working_idx >= 0; --working_idx){
    peak = find_peak(spectrum->peaks, working_idx);
    mz = get_peak_location(peak);
    
    //too small go to next peak
    if(mz < (peak_mz - sp_equalize_resolution)){
      break;
    }
    
    //beyond range, break
    if(mz > (peak_mz + sp_equalize_resolution)){
      continue;
    }

    working_intensity = get_peak_intensity(peak);
    
    //check it the intensity is larger than the peak to equalize
    if(*intensity < working_intensity){
      *intensity = working_intensity;
    }
  }
}
*/

/**
 * creates the intensity sum array with in the given resolution 
 *\returns TRUE if successfully creates the sum array, else FALSE
 */
/*
BOOLEAN_T create_intensity_sum_array(
  SCORER_T* scorer,        ///< the scorer object -in
  SPECTRUM_T* spectrum ///< the spectrum to query the intensity sum -in
  //float resolution      ///< how close are the intensities we will sum?
  )
{
  int peak_idx = 0;
  float intensity = 0;
  PEAK_T* peak = NULL;
  float peak_mz = 0;

  //check if sum array already exist
  if(spectrum->sum_array_exist){
    carp(CARP_WARNING, "intensity sum array alread existed, recreating intensity sum array");
    free(spectrum->intensity_sum_array);
  }
  
  //set sum array bounderies, and other needed variables
  float min_mz = get_scorer_sp_min_mz(scorer);
  float max_mz = get_scorer_sp_max_mz(scorer);
  float sp_array_resolution = get_scorer_sp_array_resolution(scorer);
  float sp_equalize_resolution = get_scorer_sp_equalize_resolution(scorer);

  //create the sp_array
  int array_size = (max_mz - min_mz)/sp_array_resolution;
  spectrum->intensity_sum_array = (float*)mycalloc(array_size, sizeof(float));
  spectrum->sp_sum_resolution = get_scorer_sp_sum_resolution(scorer);

  //iterate over all peaks and add the intensity to the correct array sum bins
  for(; peak_idx < spectrum->num_peaks; ++peak_idx){
    peak = find_peak(spectrum->peaks, peak_idx);
    intensity = get_peak_intensity(peak);
    peak_mz = get_peak_location(peak);
    
    //set the peak intensity to the highest of it's +/- 1u peaks
    equalize_peak(spectrum, peak_idx, &intensity, peak_mz, sp_equalize_resolution);

    //add to sum array
    add_intensity_to_sum_array(scorer, spectrum, array_size, intensity, get_peak_location(peak));
  }
  
  //now we have create the sum array
  spectrum->sum_array_exist = TRUE;
  
  return TRUE;
}
*/

/**
 * \returns The sum of intensities within 'resolution' of 'mz' in 'spectrum'
 * NOTE: Chris, this should lazily create the data structures within the
 * spectrum object that it needs.
 */
/*
float get_nearby_intensity_sum(
  SCORER_T* scorer,        ///< the scorer object -in                          
  SPECTRUM_T* spectrum, ///< the spectrum to query the intensity sum -in
  float mz             ///< the mz of the peak around which to sum intensities
  //float resolution      ///< how close are the intensities we will sum?
  )
{
  // check if sum array already exist
  if(!spectrum->sum_array_exist){
    carp(CARP_INFO, "creating intensity sum array");
    create_intensity_sum_array(scorer, spectrum);
  }
  // check if resolution matches the resolution wanted
  // if doesn't match, create a new intensity sum_array
  else if(spectrum->sp_sum_resolution != get_scorer_sp_sum_resolution(scorer)){//FIXME
    carp(CARP_INFO, "existing array doesn't support the query resolution: %.2f", get_scorer_sp_sum_resolution(scorer));
    create_intensity_sum_array(scorer, spectrum);
  }

  //get hash index into sum_array
  int hash_index = hash_sum_array(scorer, mz);

  //return the sum of the intensities
  return spectrum->intensity_sum_array[hash_index];
}
*/

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

