/*************************************************************************//**
 * \file spectrum.cpp
 * AUTHOR: Chris Park
 * CREATE DATE:  June 22 2006
 * DESCRIPTION: code to support working with spectra
 * REVISION: $Revision: 1.71 $
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
#include "carp.h"

#include "DelimitedFile.h"
#include "./MSToolkit/Spectrum.h"

//There is a name clash with MS2 in MSToolkit, so can't use the namespace decl here
//using namespace MSToolkit;


/**
 * \define constants
 */
#define MAX_PEAKS 4000
#define MZ_TO_PEAK_ARRAY_RESOLUTION 5 // i.e. 0.2 m/z unit
#define MAX_PEAK_MZ 5000
#define MAX_CHARGE 6
#define MAX_I_LINES 2 // number of 'I' lines albe to parse for one spectrum object
#define MAX_D_LINES 2 // number of 'D' lines albe to parse for one spectrum object

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
struct spectrum{
  int              first_scan;    ///< The number of the first scan
  int              last_scan;     ///< The number of the last scan
  int              id;            ///< A unique identifier
                                  // FIXME, this field is not set when parsing
  SPECTRUM_TYPE_T  spectrum_type; ///< The type of spectrum. 
  FLOAT_T            precursor_mz;  ///< The m/z of precursor (MS-MS spectra)
  int*             possible_z;    ///< The possible charge states of this spectrum
  int              num_possible_z;///< The number of possible charge states of this spectrum
  PEAK_T*          peaks;         ///< The spectrum peaks
  FLOAT_T            min_peak_mz;   ///< The minimum m/z of all peaks
  FLOAT_T            max_peak_mz;   ///< The maximum m/z of all peaks
  int              num_peaks;     ///< The number of peaks
  double           total_energy;  ///< The sum of intensities in all peaks
  char*            filename;      ///< Optional filename
  char*            i_lines[MAX_I_LINES]; ///< store i lines, upto MAX_I_LINES
  char*            d_lines[MAX_D_LINES]; ///< store d lines, upto MAX_D_LINES 
  BOOLEAN_T        has_peaks;  ///< Does the spectrum contain peak information
  BOOLEAN_T        sorted_by_mz; ///< Are the spectrum peaks sorted by m/z...
  BOOLEAN_T        sorted_by_intensity; ///< ... or by intensity?
  BOOLEAN_T        has_mz_peak_array; ///< Is the mz_peak_array populated.
  PEAK_T**         mz_peak_array;  ///< Allows rapid peak retrieval by mz.
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
 * Parses the 'D' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_D_line(
  SPECTRUM_T* spectrum, ///< place D line into this spectrum -out  
  char* line  ///< 'D' line to parse -in
);

/*
 * Parses the 'I' line of the a spectrum
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_I_line(
  SPECTRUM_T* spectrum, ///< place I line into this spectrum -out  
  char* line  ///< 'I' line to parse -in
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
  int line_idx;
  SPECTRUM_T* fresh_spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));
  fresh_spectrum->possible_z = (int*)mymalloc(sizeof(int) * MAX_CHARGE);
  fresh_spectrum->peaks = allocate_peak_array(MAX_PEAKS);
  fresh_spectrum->num_peaks = 0;
  
  // initialize D lines
  for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
    fresh_spectrum->d_lines[line_idx] = NULL;
  }

  // initialize I lines
  for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
    fresh_spectrum->i_lines[line_idx] = NULL;
  }
  
  fresh_spectrum->has_peaks = FALSE;

  return fresh_spectrum;
}

/**
 * \returns A new spectrum object, populated with the user specified parameters.
 */ 

SPECTRUM_T* new_spectrum(
  int               first_scan,         ///< The number of the first scan -in
  int               last_scan,          ///< The number of the last scan -in
  SPECTRUM_TYPE_T   spectrum_type,      ///< The type of spectrum. -in
  FLOAT_T             precursor_mz,       ///< The m/z of the precursor (for MS-MS spectra) -in
  int*              possible_z,         ///< The possible charge states of this spectrum  -in
  int               num_possible_z,     ///< The number of possible charge states of this spectrum  -in
  char*             filename)           ///< Optional filename -in
{
  SPECTRUM_T* fresh_spectrum = allocate_spectrum();
  fresh_spectrum->first_scan = first_scan;
  fresh_spectrum->last_scan = last_scan;
  fresh_spectrum->spectrum_type = spectrum_type;
  fresh_spectrum->precursor_mz = precursor_mz;
  fresh_spectrum->sorted_by_mz = FALSE;
  fresh_spectrum->has_mz_peak_array = FALSE;
  fresh_spectrum->sorted_by_intensity = FALSE;
  fresh_spectrum->peaks = NULL;
  fresh_spectrum->num_peaks = 0;
  fresh_spectrum->mz_peak_array = NULL;
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
  int line_idx;
  
  // only non post_process spectrum has these features to free
  if(spectrum->has_peaks){
    free(spectrum->possible_z);
    free(spectrum->peaks);
    free(spectrum->filename);
    
    // free D lines
    for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
      if(spectrum->d_lines[line_idx] != NULL){
        free(spectrum->d_lines[line_idx]);
      }
      else{
        break;
      }
    }
    
    // free I lines
    for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
      if(spectrum->i_lines[line_idx] != NULL){
        free(spectrum->i_lines[line_idx]);
      }
      else{
        break;
      }
    }    
  }
  
  if(spectrum->has_mz_peak_array){
    free(spectrum->mz_peak_array);
  }
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
  int num_d_index = 0;
  int num_i_index = 0;
  int num_peak_index = 0;

  fprintf(file, "S\t%06d\t%06d\t%.2f\n", 
         spectrum->first_scan,
         spectrum->last_scan,
         spectrum->precursor_mz);

  // print 'I' line
  for(; num_i_index < MAX_I_LINES; ++num_i_index){
    if(spectrum->i_lines[num_i_index] == NULL){
      break;
    }

    fprintf(file, "%s", spectrum->i_lines[num_i_index]);
  }
  
  // print 'Z', 'D' line
  for(; num_z_index < spectrum->num_possible_z; ++num_z_index){

    // print 'Z' line
    fprintf(file, "Z\t%d\t%.2f\n", spectrum->possible_z[num_z_index],
            get_spectrum_singly_charged_mass(spectrum,
                                             spectrum->possible_z[num_z_index]));

    // are there any 'D' lines to print?
    if(num_d_index < MAX_D_LINES){
      if(spectrum->d_lines[num_d_index] != NULL){
        fprintf(file, "%s", spectrum->d_lines[num_d_index]);
      }
      ++num_d_index;
    }
  }
  
  // print peaks
  for(; num_peak_index < spectrum->num_peaks; ++num_peak_index){
    fprintf(file, "%.2f %.13f\n", 
            get_peak_location(find_peak(spectrum->peaks, num_peak_index)),
            get_peak_intensity(find_peak(spectrum->peaks, num_peak_index)));
  }
}

/**
 * Prints a spectrum with the given intensities instead of the
 * observed peaks.  Assumes intensities are in m/z bins from 0 to
 * max_mz_bin.  Only prints non-zero intensities.
 */
void print_spectrum_processed_peaks(
  SPECTRUM_T* spectrum, ///< the spectrum to print 
  int charge,       ///< print at this charge state
  FLOAT_T* intensities, ///< intensities of new peaks
  int max_mz_bin,       ///< num_bins in intensities
  FILE* file){          ///< print to this file

  int i_idx = 0;
  int z_idx = 0;
  int d_idx = 0;

  // print S line
  fprintf(file, "S\t%06d\t%06d\t%.2f\n", 
         spectrum->first_scan,
         spectrum->last_scan,
         spectrum->precursor_mz);

  // print I line(s)
  for(i_idx = 0; i_idx < MAX_I_LINES; i_idx++){
    if(spectrum->i_lines[i_idx] == NULL){
      break;
    }
    fprintf(file, "%s", spectrum->i_lines[i_idx]);
  }

  // print 'Z', 'D' line
  for(z_idx = 0; z_idx < spectrum->num_possible_z; ++z_idx){

    // print 'Z' line
    if( charge == 0  || charge == spectrum->possible_z[z_idx] ){
      fprintf(file, "Z\t%d\t%.2f\n", spectrum->possible_z[z_idx],
              get_spectrum_singly_charged_mass(spectrum,
                                               spectrum->possible_z[z_idx]));
    }
    // are there any 'D' lines to print?
    if(d_idx < MAX_D_LINES){
      if(spectrum->d_lines[d_idx] != NULL){
        fprintf(file, "%s", spectrum->d_lines[d_idx]);
      }
      ++d_idx;
    }
  }

  // print peaks
  int bin_idx = 0;
  for(bin_idx = 0; bin_idx < max_mz_bin; bin_idx++){
    if( intensities[bin_idx] != 0 ){
      fprintf(file, "%d %.4f\n", bin_idx, intensities[bin_idx]); 
    }
  }
  return;
}

/**
 * Prints a spectrum object to file in sqt format.
 */
void print_spectrum_sqt(
  SPECTRUM_T* spectrum, ///< spectrum to print -in
  FILE* file,           ///< output file to print at -out
  int num_matches,      ///< number of peptides compared to this spec -in
  int charge            ///< charge used for the search -in
  ){

  int precision = get_int_parameter("precision");
  char format[64];
  sprintf(format, 
          "S\t%%d\t%%d\t%%d\t%%.%if\t%%s\t%%.%if\t%%.%if\t%%.%if\t%%d\n", 
          precision, precision, precision, precision);
  //<first scan><last scan><charge><precursor m/z><# sequence match>
  //fprintf(file, "S\t%d\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%d\n", 
  fprintf(file, format,
          get_spectrum_first_scan(spectrum), 
          get_spectrum_last_scan(spectrum),
          charge, 
          0.0, // FIXME dummy <process time>
          "server", // FIXME dummy <server>
          // get_spectrum_precursor_mz(spectrum), // TODO this should be M+H+
          get_spectrum_neutral_mass(spectrum, charge), //this is used in search
          0.0, // FIXME dummy <total intensity>
          0.0, // FIXME dummy <lowest sp>
          num_matches);
  
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
  int line_idx;

  // copy each varible
  set_spectrum_first_scan(dest,get_spectrum_first_scan(src));
  set_spectrum_last_scan(dest,get_spectrum_last_scan(src));
  set_spectrum_id(dest,get_spectrum_id(src));
  set_spectrum_spectrum_type(dest,get_spectrum_spectrum_type(src));
  set_spectrum_precursor_mz(dest,get_spectrum_precursor_mz(src));
  
  // copy possible_z
  possible_z = get_spectrum_possible_z(src);
  set_spectrum_possible_z(dest,possible_z, 
                          get_spectrum_num_possible_z(src));
  free(possible_z);
  
  // copy filename
  new_filename = get_spectrum_filename(src);
  set_spectrum_filename(dest, new_filename);
  free(new_filename);
  
  // copy 'D', 'I' lines
  for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
    if(src->d_lines[line_idx] != NULL){
      dest->d_lines[line_idx] = my_copy_string(src->d_lines[line_idx]);
    }
    else{
      break;
    }
  }
  
  // copy 'D', 'I' lines
  for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
    if(src->i_lines[line_idx] != NULL){
      dest->i_lines[line_idx] = my_copy_string(src->i_lines[line_idx]);
    }
    else{
      break;
    }
  }

  // copy each peak
  for(; num_peak_index < get_spectrum_num_peaks(src); ++num_peak_index){
    add_peak_to_spectrum(dest, get_peak_intensity(find_peak(src->peaks, num_peak_index)),
                         get_peak_location(find_peak(src->peaks, num_peak_index))); 
  }
}


/**
 * Parses a spectrum from an .mgf file
 * \returns TRUE if success. FALSE is failure.
 * 'I'
 * Skips Header line "H"
 * FIXME if need to read 'H', header line, does not parse ID
 */
BOOLEAN_T parse_spectrum_file_mgf(
  SPECTRUM_T* spectrum, ///< spectrum to parse the information into -out
  FILE* file, ///< the input file stream -in
  char* filename ///< filename of the spectrum, should not free -in
  )
{
  //long file_index = ftell(file); // stores the location of the current working line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FLOAT_T location_mz;
  FLOAT_T intensity;

  BOOLEAN_T begin_found = FALSE;
  BOOLEAN_T title_found = FALSE;
  BOOLEAN_T charge_found = FALSE;
  BOOLEAN_T pepmass_found = FALSE;
  BOOLEAN_T peaks_found = FALSE;
  BOOLEAN_T end_found = FALSE;

  carp(CARP_DEBUG, "parsing MGF Scan");

  while( (line_length = getline(&new_line, &buf_length, file)) != -1){
    //scan until BEGIN IONS
    if (strncmp(new_line, "BEGIN IONS", 10) == 0) {
      begin_found = TRUE;
      break;
    }
  }

  if (!begin_found) {
    carp(CARP_DEBUG,"Couldn't find any more scans");
    return FALSE;
  }
  //scan for the header fields
  while( (line_length = getline(&new_line, &buf_length, file)) != -1){
    if (strncmp(new_line, "TITLE=",6) == 0) {
      title_found = TRUE;
      int first_scan;
      int last_scan;
      //parse the title line
      //assume the format of title line is title.first_scan.last_scan.charge.dta
      char* period_title_ptr = index(new_line,'.');
      char* period_first_scan_ptr = index(period_title_ptr+1,'.');
      char* period_last_scan_ptr = index(period_first_scan_ptr+1,'.');
      //char* period_charge_ptr = index(period_last_scan_ptr+1,'.');

      *period_title_ptr = '\0';
      carp(CARP_DETAILED_DEBUG, "title:%s", new_line);
      
      *period_first_scan_ptr = '\0';
      carp(CARP_DETAILED_DEBUG, "first scan:%s", (period_title_ptr+1));
      first_scan = atoi((period_title_ptr+1));

      *period_last_scan_ptr = '\0';
      carp(CARP_DETAILED_DEBUG, "last scan:%s", (period_first_scan_ptr+1));
      last_scan = atoi((period_first_scan_ptr+1));

      carp(CARP_DETAILED_DEBUG, "first_scan:%d", first_scan);
      carp(CARP_DETAILED_DEBUG, "last_scan:%d", last_scan);

      set_spectrum_first_scan(spectrum, first_scan);
      set_spectrum_last_scan(spectrum, last_scan);
      set_spectrum_spectrum_type(spectrum, MS2);
    } else if (strncmp(new_line, "CHARGE=",7) == 0) {
      //parse the charge line
      int charge;
      char* plus_index = index(new_line,'+');
      *plus_index = '\0';
      carp(CARP_DETAILED_DEBUG,"Parsing %s",(new_line+7));
      charge = atoi(new_line+7);

      carp(CARP_DETAILED_DEBUG, "charge:%d", charge);

      add_possible_z(spectrum, charge);

      charge_found = TRUE;
    } else if (strncmp(new_line, "PEPMASS=",8) == 0) {
      //parse the pepmass line
      FLOAT_T pepmass;
      carp(CARP_DETAILED_DEBUG, "Parsing %s",(new_line+8));
      pepmass = atof(new_line+8);
      carp(CARP_DETAILED_DEBUG, "pepmass:%f",pepmass);
      //TODO - check to see if this is correct.
      set_spectrum_precursor_mz(spectrum, pepmass);
      pepmass_found = TRUE;
    } else if (isdigit(new_line[0])) {
      //no more header lines, peak information is up
      peaks_found = TRUE;
      break;
    } else if (strcmp(new_line, "END IONS")) {
      //we found the end of the ions without any peaks.
      carp(CARP_WARNING,"No peaks found for mgf spectrum");
      return TRUE;
    }
  }

  //TODO check to make sure we gleaned the information from
  //the headers.

  //parse peak information
  do {
    if (strncmp(new_line, "END IONS", 8) == 0) {
      //we are done parsing this charged spectrum.
      end_found = TRUE;
      break;
    }
    #ifdef USE_DOUBLES
    else if(sscanf(new_line,"%lf %lf", &location_mz, &intensity) == 2){
    #else
    else if(sscanf(new_line,"%f %f", &location_mz, &intensity) == 2){
    #endif
      carp(CARP_DETAILED_DEBUG,"adding peak %f %f",location_mz, intensity);
      //add the peak to the spectrum object
      add_peak_to_spectrum(spectrum, intensity, location_mz);
    } else {
      //file format error.
      carp(CARP_ERROR,
      "File format error\n"
      "At line: %s",
      new_line);
    }
  } while( (line_length = getline(&new_line, &buf_length, file)) != -1);

  if (end_found) {
    //we successfully parsed this spectrum.
    spectrum -> has_peaks = TRUE;
    set_spectrum_new_filename(spectrum, filename);  
    return TRUE;
  } else {
    //something happened, bomb.
    return FALSE;
  }
}
 

/**
 * Parses a spectrum from file.
 * \returns TRUE if success. FALSE is failure.
 * 'I'
 * Skips Header line "H"
 * FIXME if need to read 'H', header line, does not parse ID
 */
BOOLEAN_T parse_spectrum_file(
  SPECTRUM_T* spectrum, ///< spectrum to parse the information into -out
  FILE* file, ///< the input file stream -in
  char* filename ///< filename of the spectrum, should not free -in
  )
{

  if (get_boolean_parameter("use-mgf")) {
    return parse_spectrum_file_mgf(spectrum, file, filename);
  }

  long file_index = ftell(file); // stores the location of the current working line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  FLOAT_T location_mz;
  FLOAT_T intensity;
  BOOLEAN_T record_S = FALSE; // check's if it read S line
  BOOLEAN_T record_Z = FALSE; // check's if it read Z line
  BOOLEAN_T start_add_peaks = FALSE; // check's if it started reading peaks
  BOOLEAN_T file_format = FALSE; // is the file format correct so far
  
  FLOAT_T test_float;
  char test_char;
  
  while( (line_length = getline(&new_line, &buf_length, file)) != -1){
    // skip header line
    // if(new_line[0] == 'H'){
    //  file_index = ftell(file);
    //  continue;
    // }
    // checks if 'S' is not the first line
    if((!record_S || (record_S && start_add_peaks)) && 
            (new_line[0] == 'Z' ||  
             new_line[0] == 'I' ||
             new_line[0] == 'D' )){
      file_format = FALSE;
      carp(
        CARP_ERROR, 
        "Incorrect order of line (S,Z, Peaks)\n"
        "At line: %s", 
        new_line
      );
      break; // File format incorrect
    }
    // Reads the 'S' line
    else if(new_line[0] == 'S' && !record_S){
      record_S = TRUE;
      if(!parse_S_line(spectrum, new_line, buf_length)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    // Reads the 'Z' line 
    else if(new_line[0] == 'Z'){
      record_Z = TRUE;
      if(!parse_Z_line(spectrum, new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }

    // Reads the 'D' line 
    else if(new_line[0] == 'D'){
      if(!parse_D_line(spectrum, new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }

    // Reads the 'I' line 
    else if(new_line[0] == 'I'){
      if(!parse_I_line(spectrum, new_line)){
        file_format = FALSE;
        break; // File format incorrect
      }
    }
    
    // Stops, when encounters the start of next spectrum 'S' line
    else if(new_line[0] == 'S' && start_add_peaks){ // start of next spectrum
      break;
    }

    // *****parse peak line******
    else if(new_line[0] != 'Z' &&  
            new_line[0] != 'I' &&
            new_line[0] != 'D' &&
            new_line[0] != '\n')
      {
        // checks if the peaks are in correct order of lines
        if((!record_Z || !record_S)){
          file_format = FALSE;
          carp(
            CARP_ERROR,
            "Incorrect order of line (S,Z, Peaks)\n"
            "At line: %s", 
            new_line
          );
          break; // File format incorrect
        }
        // check for peak line format
        #ifdef USE_DOUBLES
        else if((sscanf(new_line,"%lf %lf %lf",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_float) > 2)||
                (sscanf(new_line,"%lf %lf %c",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_char) > 2)||
                (sscanf(new_line,"%lf %lf",// test format:peak line has less than 2 fields
                        &test_float, &test_float) != 2)){
        #else
        else if((sscanf(new_line,"%f %f %f",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_float) > 2)||
                (sscanf(new_line,"%f %f %c",// test format:peak line has more than 2 fields
                        &test_float, &test_float, &test_char) > 2)||
                (sscanf(new_line,"%f %f",// test format:peak line has less than 2 fields
                        &test_float, &test_float) != 2)){
        #endif
          file_format = FALSE;
          carp(
            CARP_ERROR,
            "Incorrect peak line\n"
            "At line: %s", 
            new_line
          );
          break; // File format incorrect
        }
        // Reads the 'peak' lines, only if 'Z','S' line has been read
        #ifdef USE_DOUBLES
        else if(record_Z && record_S &&
                (sscanf(new_line,"%lf %lf", &location_mz, &intensity) == 2)){
        #else
        else if(record_Z && record_S &&
                (sscanf(new_line,"%f %f", &location_mz, &intensity) == 2)){
        #endif
          file_format = TRUE;
          start_add_peaks = TRUE;
          add_peak_to_spectrum(spectrum, intensity, location_mz);
        }
      }
    // *************************
    file_index = ftell(file); // updates the current working line location
  }

  // now we have peak information
  spectrum->has_peaks = TRUE;
    
  // set the file pointer back to the start of the next 's' line
  fseek(file, file_index, SEEK_SET);
  myfree(new_line);
  
  // set filename of empty spectrum
  set_spectrum_new_filename(spectrum, filename);
  
  // No more spectrum in .ms file
  if(!record_S && !file_format){
    return FALSE;
  }
  
  // File format incorrect
  if(!file_format){ 
    carp(CARP_ERROR, "incorrect file format\n");
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
   FLOAT_T precursor_mz;
   FLOAT_T test_float;
   char test_char;

   // deletes empty space & 0
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
   // deletes empty space & zeros
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
   #ifdef USE_DOUBLES
   if ( (sscanf(spliced_line,"%lf %lf %lf %lf",// test format:S line has more than 3 fields
                &test_float, &test_float, &test_float, &test_float) > 3) ||
        (sscanf(spliced_line,"%lf %lf %lf %c",// test format:S line has more than 3 fields 
                &test_float, &test_float, &test_float, &test_char) > 3) ||
        (sscanf(spliced_line,"%i %i %lf", // S line is parsed here
               &first_scan, &last_scan, &precursor_mz) != 3)) {
   #else
   if ( (sscanf(spliced_line,"%f %f %f %f",// test format:S line has more than 3 fields
                &test_float, &test_float, &test_float, &test_float) > 3) ||
        (sscanf(spliced_line,"%f %f %f %c",// test format:S line has more than 3 fields 
                &test_float, &test_float, &test_float, &test_char) > 3) ||
        (sscanf(spliced_line,"%i %i %f", // S line is parsed here
               &first_scan, &last_scan, &precursor_mz) != 3)) {
   #endif
     carp(CARP_ERROR,"Failed to parse 'S' line:\n %s",line);
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
   FLOAT_T m_h_plus;
   FLOAT_T test_float;
   char test_char;

   // check if Z line is in correct format
   #ifdef USE_DOUBLES
   if( ((tokens =  // test format: Z line has less than 3 fields
         sscanf(line, "%c %lf %lf", &test_char, &test_float, &test_float)) < 3) ||
       ((tokens =   // test format: Z line has more than 3 fields
         sscanf(line, "%c %lf %lf %lf", &test_char, &test_float, &test_float, &test_float)) >  3) ||
       ((tokens =  // test format: Z line has more than 3 fields
         sscanf(line, "%c %lf %lf %c", &test_char, &test_float, &test_float, &test_char)) >  3) ||
       (tokens = // Z line is parsed here
        sscanf(line, "%c %d %lf", &line_name, &charge, &m_h_plus)) != 3){
   #else
   if( ((tokens =  // test format: Z line has less than 3 fields
         sscanf(line, "%c %f %f", &test_char, &test_float, &test_float)) < 3) ||
       ((tokens =   // test format: Z line has more than 3 fields
         sscanf(line, "%c %f %f %f", &test_char, &test_float, &test_float, &test_float)) >  3) ||
       ((tokens =  // test format: Z line has more than 3 fields
         sscanf(line, "%c %f %f %c", &test_char, &test_float, &test_float, &test_char)) >  3) ||
       (tokens = // Z line is parsed here
        sscanf(line, "%c %d %f", &line_name, &charge, &m_h_plus)) != 3){
   #endif
     carp(CARP_ERROR,"Failed to parse 'Z' line:\n %s",line);
     return FALSE;
   }  

   return add_possible_z(spectrum, charge);
 }

 /**
  * Adds a possible charge(z) to the spectrum.
  * Must not exceed the MAX_CARGE capacity
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
  * FIXME currently does not parse D line, just copies the entire line
  * Parses the 'D' line of the a spectrum
  * \returns TRUE if success. FALSE is failure.
  */
 BOOLEAN_T parse_D_line(
   SPECTRUM_T* spectrum, ///< place D line into this spectrum -out  
   char* line  ///< 'D' line to parse -in
   )
 {
   int line_idx;
   int length = strlen(line)+1;
   char* d_line = (char*)mycalloc(length, sizeof(char));

   strncpy(d_line, line, length-3);
   d_line[length-2] = '\0';
   d_line[length-3] = '\n';
   
   // find empty spot D lines
   for(line_idx = 0; line_idx < MAX_D_LINES; ++line_idx){
     // check for empty space
     if(spectrum->d_lines[line_idx] == NULL){
       spectrum->d_lines[line_idx] = d_line;
       break;
     }
   }

   // check if added new d line to spectrum
   if(line_idx == MAX_D_LINES){
     free(d_line);
     carp(CARP_WARNING, "no more space for additional D lines, max: %d", MAX_D_LINES);
   }

   return TRUE;
 }

 /**
  * FIXME currently does not parse I line, just copies the entire line
  * Parses the 'I' line of the a spectrum
  * \returns TRUE if success. FALSE is failure.
  */
 BOOLEAN_T parse_I_line(
   SPECTRUM_T* spectrum, ///< place I line into this spectrum -out  
   char* line  ///< 'I' line to parse -in
   )
 {
   int line_idx;
   int length = strlen(line)+1;
   char* i_line = (char*)mycalloc(length, sizeof(char));

   strncpy(i_line, line, length-3);
   i_line[length-2] = '\0';
   i_line[length-3] = '\n';

  // find empty spot I lines
  for(line_idx = 0; line_idx < MAX_I_LINES; ++line_idx){
    // check for empty space
    if(spectrum->i_lines[line_idx] == NULL){
      spectrum->i_lines[line_idx] = i_line;
      break;
    }
  }

  // check if added new i line to spectrum
  if(line_idx == MAX_I_LINES){
    free(i_line);
    carp(CARP_WARNING, "no more space for additional I lines, max: %d", MAX_I_LINES);
  }

  return TRUE;
}

/**
  * Parses a spectrum from a MSToolkit::Spectrum
  * \returns TRUE if success. FALSE is failure.
  */
BOOLEAN_T parse_spectrum_spectrum(
  SPECTRUM_T* spectrum, ///< spectrum to parse the information into -out
  void* mst_spectrum, ///< the input MSToolkit spectrum -in
  char* filename ///< filename of the spectrum, should not free -in
  ) {


  MSToolkit::Spectrum* mst_real_spectrum = (MSToolkit::Spectrum*)mst_spectrum;

  //set first_scan, last_scan, and precursor_mz.
  set_spectrum_first_scan(spectrum, mst_real_spectrum -> getScanNumber());
  set_spectrum_last_scan(spectrum, mst_real_spectrum -> getScanNumber());
  set_spectrum_precursor_mz(spectrum, mst_real_spectrum -> getMZ());

  // set filename of empty spectrum
  set_spectrum_new_filename(spectrum, filename);

  spectrum -> spectrum_type = MS2; //assume MS2
  

  //add all peaks.
  
  int peak_idx;
  for (peak_idx=0;peak_idx < mst_real_spectrum -> size();peak_idx++) {
    add_peak_to_spectrum(spectrum,
      mst_real_spectrum -> at(peak_idx).intensity,
      mst_real_spectrum -> at(peak_idx).mz);
  }
  
  //add possible charge states.
  if(  mst_real_spectrum->sizeZ() > 0 ){
    for (int z_idx = 0; z_idx < mst_real_spectrum -> sizeZ(); z_idx++) {
      add_possible_z(spectrum, mst_real_spectrum -> atZ(z_idx).z);
    }
  } else { // if no charge states detected, decide based on spectrum
    int charge = choose_charge(spectrum->precursor_mz,
                               spectrum->peaks,
                               spectrum->num_peaks);

    // add either +1 or +2, +3
    if( charge == 1 ){
      add_possible_z(spectrum, 1);
    } else if( charge == 0 ){
      add_possible_z(spectrum, 2);
      add_possible_z(spectrum, 3);
    } else {
      carp(CARP_ERROR, "Could not determine charge state for spectrum %d.", 
           spectrum->first_scan);
    }
  }

  spectrum -> has_peaks = TRUE;

  return TRUE;
}

/**
 * Adds a peak to the spectrum given a intensity and location
 * calls update_spectrum_fields to update num_peaks, min_peak ...
 */
BOOLEAN_T add_peak_to_spectrum(
  SPECTRUM_T* spectrum,///< spectrum to add the peak to -out 
  FLOAT_T intensity, ///< the intensity of peak to add -in
  FLOAT_T location_mz ///< the location of peak to add -in
  )
{
  if(spectrum->num_peaks < MAX_PEAKS){  // FIXME change it to be dynamic
    set_peak_intensity(find_peak(spectrum->peaks, spectrum->num_peaks),
                       intensity);
    set_peak_location(find_peak(spectrum->peaks, spectrum->num_peaks),
                      location_mz);
    update_spectrum_fields(spectrum, intensity, location_mz);
    spectrum->has_peaks = TRUE;
    return TRUE;
  }

  return FALSE;
}

void populate_mz_peak_array(
  SPECTRUM_T* spectrum
  )
{
  if (spectrum->has_mz_peak_array == TRUE){
    return;
  }
  
  int array_length = MZ_TO_PEAK_ARRAY_RESOLUTION * MAX_PEAK_MZ;
  PEAK_T** mz_peak_array = (PEAK_T**)mymalloc(array_length * sizeof(PEAK_T*));
  int peak_idx;
  for (peak_idx = 0; peak_idx < array_length; peak_idx++){
    mz_peak_array[peak_idx] = NULL;
  }
  PEAK_ITERATOR_T* peak_iterator = new_peak_iterator(spectrum);
  PEAK_T* peak = NULL;
  while(peak_iterator_has_next(peak_iterator)){
    peak = peak_iterator_next(peak_iterator);
    FLOAT_T peak_mz = get_peak_location(peak);
    int mz_idx = (int) (peak_mz * MZ_TO_PEAK_ARRAY_RESOLUTION);
    if (mz_peak_array[mz_idx] != NULL){
      carp(CARP_INFO, "Peak collision at mz %.3f = %i", peak_mz, mz_idx);
      if(get_peak_intensity(mz_peak_array[mz_idx])< get_peak_intensity(peak)){
        mz_peak_array[mz_idx] = peak;
      }
    } else {
      mz_peak_array[mz_idx] = peak; 
    }
  }
  spectrum->mz_peak_array = mz_peak_array;
  spectrum->has_mz_peak_array = TRUE;
  free_peak_iterator(peak_iterator);
}

/**
 * \returns The closest intensity within 'max' of 'mz' in 'spectrum'
 * NULL if no peak.
 * This should lazily create the data structures within the
 * spectrum object that it needs.
 * TODO: reimplement with faster peak lookup
 */
PEAK_T* get_nearest_peak(
  SPECTRUM_T* spectrum, ///< the spectrum to query the intensity sum -in
  FLOAT_T mz, ///< the mz of the peak around which to sum intensities -in
  FLOAT_T max ///< the maximum distance to get intensity -in
  )
{
  populate_mz_peak_array(spectrum); // for rapid peak lookup by mz

  FLOAT_T min_distance = BILLION;
  int min_mz_idx = (int)((mz - max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  min_mz_idx = min_mz_idx < 0 ? 0 : min_mz_idx;
  int max_mz_idx = (int)((mz + max) * MZ_TO_PEAK_ARRAY_RESOLUTION + 0.5);
  int absolute_max_mz_idx = MAX_PEAK_MZ * MZ_TO_PEAK_ARRAY_RESOLUTION - 1;
  max_mz_idx = max_mz_idx > absolute_max_mz_idx 
    ? absolute_max_mz_idx : max_mz_idx;
  PEAK_T* peak = NULL;
  PEAK_T* nearest_peak = NULL;
  int peak_idx;
  for (peak_idx=min_mz_idx; peak_idx < max_mz_idx + 1; peak_idx++){
    if ((peak = spectrum->mz_peak_array[peak_idx]) == NULL){
      continue;
    }
    FLOAT_T peak_mz = get_peak_location(peak);
    FLOAT_T distance = fabs(mz - peak_mz);
    if (distance > max){
      continue;
    }
    if (distance < min_distance){
      nearest_peak = peak;
      min_distance = distance;
    }
  }
  return nearest_peak;
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
    carp(CARP_FATAL,"File %s could not be opened",filename);
    return (FALSE);  // exit(1);
  }
  // might check if spectrum is NULL?? Close file??
  if(parse_spectrum_file(spectrum, file, filename)){
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
  FLOAT_T intensity, ///< the intensity of the peak that has been added -in
  FLOAT_T location ///< the location of the peak that has been added -in
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
  // update total_energy
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
FLOAT_T get_spectrum_precursor_mz(
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
  FLOAT_T precursor_mz ///< the precursor_mz -in
  )
{
  spectrum->precursor_mz = precursor_mz;
}

/**
 * \returns the minimum m/z of all peaks
 */
FLOAT_T get_spectrum_min_peak_mz(
  SPECTRUM_T* spectrum ///< the spectrum to query min_peak_mz -in
  )
{
  return spectrum->min_peak_mz;
}

/**
 * \returns the maximum m/z of all peaks
 */
FLOAT_T get_spectrum_max_peak_mz(
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
  
  int filename_length = strlen(spectrum->filename) +1; // +\0
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
  int filename_length = strlen(filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  spectrum->filename =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * \returns the number of possible charge states of this spectrum
 */
int get_num_possible_z(
  SPECTRUM_T* spectrum ///< the spectrum to query possible z -in
  )
{
  return spectrum->num_possible_z;
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
 * \returns a pointer to an array of the possible charge states of this spectrum
 * User must NOT free this or alter, not a copy
 * number of possible charge states can be gained by 
 * the get_num_possible_z function.
 */
int* get_spectrum_possible_z_pointer(
  SPECTRUM_T* spectrum  ///< the spectrum to query possible z -in
  )
{
  return spectrum->possible_z;
}

/**
 *  Considers all parameters and allocates an array of 
 *  charges that should be searched.  Returns the number of charges
 *  in that array.  Protects get_spectrum_possible_z_pointer, could make it
 *  private.
 */ 
int get_charges_to_search(SPECTRUM_T* spectrum, int** select_charge_array){

  int total_charges = spectrum->num_possible_z;
  int* all_charge_array = spectrum->possible_z;

  int param_charge = 0;
  const char* charge_str = get_string_parameter_pointer("spectrum-charge");
  int i=0;

  // Return full array of charges
  if( strcmp( charge_str, "all") == 0){

    *select_charge_array = (int*)mymalloc(sizeof(int) * total_charges);
    for(i=0; i < total_charges; i++){
      (*select_charge_array)[i] = all_charge_array[i];
    }
    return total_charges;
  }
  // else return one charge

  param_charge = atoi(charge_str);

  if( (param_charge < 1) || (param_charge > MAX_CHARGE) ){
    carp(CARP_FATAL, "spectrum-charge option must be 1,2,3,.. %d or 'all'.  " \
         "%s is not valid", MAX_CHARGE, charge_str);
  }


  for(i=0; i<total_charges; i++){
     
    if( all_charge_array[i] == param_charge ){ 
      *select_charge_array = (int*)mymalloc(sizeof(int));
      **select_charge_array = param_charge;
      return 1;
    }
  }

  // Else none to be searched
  *select_charge_array = NULL;
  return 0;

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
FLOAT_T get_spectrum_max_peak_intensity(
  SPECTRUM_T* spectrum  ///< the spectrum to query maximum peak intensity -in
  )
{
  int num_peak_index = 0;
  FLOAT_T max_intensity = -1;

  for(; num_peak_index < get_spectrum_num_peaks(spectrum); ++num_peak_index){
    if(max_intensity <= get_peak_intensity(find_peak(spectrum->peaks, num_peak_index))){
      max_intensity = get_peak_intensity(find_peak(spectrum->peaks, num_peak_index));
    }
  }
  return max_intensity; 
}


/**
 * \returns The mass of the charged precursor ion, according to the formula 
 * mass = m/z * charge
 */
FLOAT_T get_spectrum_mass(
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
FLOAT_T get_spectrum_neutral_mass(
  SPECTRUM_T* spectrum,  ///< the spectrum to query neutral_mass -in
  int charge ///< the charge of precursor ion -in
  )
{
  return (get_spectrum_mass(spectrum, charge) - MASS_PROTON * charge); // TESTME
}

/**
 * \returns The mass of the singly charged precursor ion, according to the formula 
 * mass = m/z * charge - (mass_H * (charge - 1))
 */
FLOAT_T get_spectrum_singly_charged_mass(
  SPECTRUM_T* spectrum,  ///< the spectrum to query charged_mass -in
  int charge ///< the charge of the precursor ion -in
  )
{
  return (get_spectrum_mass(spectrum, charge) - MASS_PROTON*(charge-1));  // TESTME
}


/**
 * serialize the spectrum in binary
 * Form,
 * <int: first_scan><int: last_scan><int: id><SPECTRUM_TYPE_T: spectrum_type>
 * <float: precursor_mz><float: retention_time>
 */
void serialize_spectrum(
  SPECTRUM_T* spectrum, ///< the spectrum to serialize -in
  FILE* file ///< output stream -out
  )
{
  // serialize the spectrum struct
  fwrite(spectrum, sizeof(SPECTRUM_T), 1, file);

}

/**
 * Parse the spectrum from the tab-delimited result file
 *\returns the parsed spectrum , else returns NULL for failed parse
 */
SPECTRUM_T* parse_spectrum_tab_delimited(
  DelimitedFile& file ///< output stream -out
  ) {

  SPECTRUM_T* spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));

  spectrum -> first_scan = file.getInteger("scan");
  spectrum -> last_scan = file.getInteger("scan");
  spectrum -> spectrum_type = MS2; //assume MS2;

  spectrum -> precursor_mz = file.getFloat("spectrum precursor m/z");
  //Is it okay to assign an individual spectrum object for each charge?
  //add_possible_z(spectrum, file.getInteger("charge")); 

  /*
  TODO : Implement these in the tab delimited file?
  spectrum -> min_peak_mz = file.getFloat("spectrum min peak mz");
  spectrum -> max_peak_mz = file.getFloat("spectrum max peak mz");
  spectrum -> num_peaks = file.getInteger("spectrum num peaks");
  spectrum -> total_energy = file.getInteger("spectrum total energy");
  */

  spectrum -> has_peaks = FALSE;
  return spectrum;

}


/**
 * Parse the spectrum from the serialized spectrum
 *\returns the parsed spectrum , else returns NULL for failed parse
 */
SPECTRUM_T* parse_spectrum_binary(
  FILE* file ///< output stream -out
  )
{
  SPECTRUM_T* spectrum = (SPECTRUM_T*)mycalloc(1, sizeof(SPECTRUM_T));
  
  // get spectrum struct
  if(fread(spectrum, (sizeof(SPECTRUM_T)), 1, file) != 1){
    carp(CARP_ERROR, "serialized file corrupted, incorrect spectrum format");
    free(spectrum);
    return NULL;
  }
  
  spectrum->has_peaks = FALSE;
  
  return spectrum;
}

/***********************************************************************
 * Normalize peak intensities so that they sum to unity.
 ***********************************************************************/
void sum_normalize_spectrum(
  SPECTRUM_T* spectrum
  )
{
  PEAK_T* peak = NULL;
  PEAK_ITERATOR_T* peak_iterator = new_peak_iterator(spectrum);
  while(peak_iterator_has_next(peak_iterator)){
    peak = peak_iterator_next(peak_iterator);
    FLOAT_T new_intensity = get_peak_intensity(peak) / spectrum->total_energy;
    set_peak_intensity(peak, new_intensity);
  }
  free_peak_iterator(peak_iterator);
}

/***********************************************************************
 * Populate peaks with rank information.
 ***********************************************************************/
void spectrum_rank_peaks(
  SPECTRUM_T* spectrum
  )
{
  PEAK_T* peak = NULL;
  PEAK_ITERATOR_T* peak_iterator = new_peak_iterator(spectrum);
  sort_peaks(spectrum->peaks, spectrum->num_peaks, _PEAK_INTENSITY);
  spectrum->sorted_by_intensity = TRUE;
  spectrum->sorted_by_mz = FALSE;
  int rank = spectrum->num_peaks;
  while(peak_iterator_has_next(peak_iterator)){
    peak = peak_iterator_next(peak_iterator);
    FLOAT_T new_rank = rank/(float)spectrum->num_peaks;
    rank--;
    set_peak_intensity_rank(peak, new_rank); 
  }
  free_peak_iterator(peak_iterator);
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
  // now we have peak information
  if(!spectrum->has_peaks){
    carp(CARP_FATAL, "Spectrum does not contain peak information");
  }

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

/**
 *  Resets the iterator to the first element
 */
void peak_iterator_reset(
  PEAK_ITERATOR_T* peak_iterator  ///< the interator for the peaks -in
  )
{
  peak_iterator->peak_index = 0;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

