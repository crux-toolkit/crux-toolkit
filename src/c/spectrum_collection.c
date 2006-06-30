/*****************************************************************************
 * \file spectrum_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * DESCRIPTION: code to support working with collection of multiple spectra
 * REVISION: $Revision: 1.7 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "spectrum.h"
#include "spectrum_collection.h" 
#include "peak.h"
#include "utils.h"
#include "unistd.h"

#define MAX_SPECTRA 40000



long binary_search_spectrum(FILE* file, int first_scan);

int match_first_scan_line(
  char* line, 
  int buf_length, 
  int query_first_scan
  );

/**
 * \struct spectrum_collection 
 */
struct spectrum_collection {
  SPECTRUM_T* spectra[MAX_SPECTRA];  ///< The spectrum peaks
  int  num_spectra;     ///< The number of spectra
  char*   filename;     ///< Optional filename
  char* comment;        ///< The spectrum_collection header lines
  BOOLEAN_T is_parsed; ///< Have we parsed all the spectra from the file?
};    

/**
 * \struct spectrum_iterator
 */
struct spectrum_iterator {
  SPECTRUM_T* spectrum; ///< The spectrum whose peaks to iterate over. 
  int  peak_index;        ///< The index of the current peak
};


/**
 * \returns An (empty) spectrum_collection object.
 */
SPECTRUM_COLLECTION_T * allocate_spectrum_collection(void){
  SPECTRUM_COLLECTION_T* collection =
    (SPECTRUM_COLLECTION_T*)mycalloc(1,sizeof(SPECTRUM_COLLECTION_T));
  collection->is_parsed = FALSE;
  return collection;
}

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does *NOT* parse all spectra in the file. 
 * This will be done lazily depending on the subsequent method
 * calls (parse_spectrum_collection get_spectrum_collection_spectrum).
 * \returns  SPECTRUM_COLLECTION_T
 */
SPECTRUM_COLLECTION_T* new_spectrum_collection(char* filename)///< The spectrum collection filename.
{
  SPECTRUM_COLLECTION_T* spectrum_collection =  allocate_spectrum_collection();

  if(access(filename, F_OK)){
    fprintf(stderr,"File %s could not be opened\n",filename);
    free(spectrum_collection);
    exit(1);
  } //FIXME check if file is empty
  
  set_spectrum_collection_filename(spectrum_collection, filename);
  
  return spectrum_collection;
}

/**
 * Frees an allocated spectrum_collection object.
 */
void free_spectrum_collection(SPECTRUM_COLLECTION_T* spectrum_collection){
  int spectrum_index = 0;
  for(;spectrum_index < spectrum_collection->num_spectra; ++spectrum_index){
    free_spectrum(spectrum_collection->spectra[spectrum_index]);
  }
  free(spectrum_collection->filename);
  free(spectrum_collection->comment);
  free(spectrum_collection);
}


/**  //////TESTME////
 * \sets the filename of the ms2 file the spectra was parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_collection_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_collection_filename(SPECTRUM_COLLECTION_T* spectrum_collection, char* filename){
  int filename_length = strlen(filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  spectrum_collection->filename =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan. Use binary search
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
BOOLEAN_T get_spectrum_collection_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< The spectrum collection
  int first_scan,      ///< The first scan of the spectrum to retrieve 
  SPECTRUM_T* spectrum ///< The (empty) allocated SPECTRUM_T object
  )
{
  FILE* file;
  long target_index;
  //check if file is still avaliable
  if ((file = fopen(spectrum_collection->filename,"r")) == NULL) {
    fprintf(stderr,"File %s could not be opened",spectrum_collection->filename);
    return (FALSE);
  }
  
  target_index = binary_search_spectrum(file, first_scan);
  // first_scan not found
  if(target_index == -1){
    fclose(file);
    return FALSE;
  }
  fseek(file, target_index, SEEK_SET);
  // failed to parse spectrum
  if(!parse_spectrum_file(spectrum, file)){
    fclose(file);
    return FALSE;
  }
  fclose(file);
  return TRUE;
}

/**
 * 
 *\returns the file position of the query spectrum begins in file
 * represented as an long integer, the return value of ftell()
 *\returns -1 if failed to find the query spectrum
 */
long binary_search_spectrum(FILE* file, int first_scan){   
  long low_index  = ftell(file); 
  long high_index;
  long mid_index;
  long working_index;
  long end_of_file_index;

  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;

  //set initial high and low position
  if(fseek(file, 1,SEEK_END) != -1){
    high_index = ftell(file);
    end_of_file_index = high_index;
  }
  else{
    fprintf(stderr, "error: file corrupted");
    return -1;
  }  
  while(low_index <= high_index){

    mid_index = (low_index + high_index)/2;
    
    if(fseek(file, mid_index,SEEK_SET) != -1){
      
      working_index = ftell(file);
      //check each line until reach 'S' line
      while((line_length =  getline(&new_line, &buf_length, file)) != -1){

        if(new_line[0] == 'S'){
          int compare = 
            match_first_scan_line(new_line, buf_length, first_scan);
          if(compare == 0){
            free(new_line);
            return working_index; //found the query match
          }
          else if(compare == 1){
            high_index = mid_index-1;
            break;
          }
          else{ //compare == -1
            low_index = mid_index+1;
            break;
          }
        }
        // store the next working line
        working_index = ftell(file);
        
        //if the High is at EOF or working_index starts looking bellow high
        if(working_index > high_index || 
           working_index == end_of_file_index-1 ){ 
          high_index = mid_index-1;
          break;
        } 
      }
    }
    else{
      fprintf(stderr, "error: file corrupted");
      return -1;
    }
  }
  free(new_line);
  return -1; //failed to find the query spectrum
}


/**
 * Parses the 'S' line of the a spectrum,
 * check if the first_scan match the query scan number
 * \returns 0 if match, 
 * \returns 1 if query_first_scan < first_scan
 * \returns -1 if query_first_scan > first_scan 
*/
int match_first_scan_line(
  char* line, 
  int buf_length, 
  int query_first_scan
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
    fprintf(stderr,"Incorrect file format\n");
    exit(1); //FIXME check if this is corrrect 
  }
  if(first_scan == query_first_scan){
    return 0;
  }
  else if(first_scan >  query_first_scan){
    return 1;
  }
  return -1; //first_scan <  query_first_scan
}




/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
