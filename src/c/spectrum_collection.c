/*************************************************************************//**
 * \file spectrum_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * DESCRIPTION: code to support working with collection of multiple spectra
 * REVISION: $Revision: 1.43 $
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include "objects.h"
#include "spectrum.h"
#include "spectrum_collection.h" 
#include "protein_index.h" 
#include "peak.h"
#include "utils.h"
#include "unistd.h"
#include "parameter.h"

#include "CMSReader.h"
#include "CSpectrum.h"

#define MAX_SPECTRA 40000 ///< max number of spectrums
#define MAX_COMMENT 1000 ///< max length of comment


/* Private functions */
long binary_search_spectrum(FILE* file, int first_scan);

int match_first_scan_line(
  char* line, 
  int buf_length, 
  int query_first_scan
  );

void queue_next_spectrum(FILTERED_SPECTRUM_CHARGE_ITERATOR_T* it);

/**
 * \struct spectrum_collection 
 * \brief A object to group together one or more spectrum objects.
 */
struct spectrum_collection {
  SPECTRUM_T* spectra[MAX_SPECTRA];  ///< The spectrum peaks
  int  num_spectra;     ///< The number of spectra
  int  num_charged_spectra; ///< The number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra)
  char* filename;     ///< Optional filename
  char comment[MAX_COMMENT];    ///< The spectrum_collection header lines
  BOOLEAN_T is_parsed; ///< Have we parsed all the spectra from the file?
};    

/**
 * \struct spectrum_iterator
 * \brief An object to iterate over the spectra within a spectrum_collection.
 */
struct spectrum_iterator {
  SPECTRUM_COLLECTION_T* spectrum_collection;///< The collection whose spectrum to iterate over. 
  int  spectrum_index; ///< The index of the current spectrum;
};

/**
 * \struct filtered_spectrum_charge_iterator
 * \brief An object to iterate over the spectra within a spectrum_collection.
 */
struct filtered_spectrum_charge_iterator {
  SPECTRUM_COLLECTION_T* spectrum_collection;///< spectra to iterate over
  BOOLEAN_T has_next;  ///< is there a spec that passes criteria
  int  spectrum_index; ///< The index of the current spectrum
  int* charges;        ///< Array of possible charges to search
  int num_charges;     ///< how many charges does the cur spec have
  int charge_index;    ///< The index of the charge of the current spectrum
  double min_mz;       ///< return only spec above this mz
  double max_mz;      ///< return only spec below this mz
  int search_charge;   ///< which z to search, 0 for all
};


/**
 * \returns An (empty) heap allocated spectrum_collection object.
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
SPECTRUM_COLLECTION_T* new_spectrum_collection(
  const char* filename///< The spectrum collection filename. -in
  )
{
  SPECTRUM_COLLECTION_T* spectrum_collection =  allocate_spectrum_collection();
  #if DARWIN
  char path_buffer[PATH_MAX];
  char* absolute_path_file =  realpath(filename, path_buffer);
  #else
  char* absolute_path_file =  realpath(filename, NULL);
  #endif
  if (absolute_path_file == NULL){
    carp(CARP_FATAL,
	   "Error from file '%s'. (%s)",
	   filename,
	   strerror(errno)); 
  }
  
  if(access(absolute_path_file, F_OK)){
    free(spectrum_collection);
    free(absolute_path_file);
    carp(CARP_FATAL,"File %s could not be opened\n", absolute_path_file);
  } // FIXME check if file is empty
  
  set_spectrum_collection_filename(spectrum_collection, absolute_path_file);
  #ifndef DARWIN
  free(absolute_path_file);
  #endif
  
  return spectrum_collection;
}

/**
 * Frees an allocated spectrum_collection object.
 */
void free_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum collection to free - in
)
{
  int spectrum_index = 0;
  for(;spectrum_index < spectrum_collection->num_spectra; ++spectrum_index){
    free_spectrum(spectrum_collection->spectra[spectrum_index]);
  }
  free(spectrum_collection->filename);
  free(spectrum_collection);
}

/**
 * Prints a spectrum_collection object to file.
 */
void print_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< spectrum_collection to print -in 
  FILE* file ///< file for output -out
  )
{
  SPECTRUM_ITERATOR_T* spectrum_iterator;
  
  fprintf(file,"comment:\n%s", spectrum_collection->comment);
  fprintf(file,".ms2 Filename: %s\n", spectrum_collection->filename);
  
  // print each spectrum
  spectrum_iterator = new_spectrum_iterator(spectrum_collection);
  while(spectrum_iterator_has_next(spectrum_iterator)){
    print_spectrum(spectrum_iterator_next(spectrum_iterator), file);
  }
  free_spectrum_iterator(spectrum_iterator);
}

/**
 * Copies spectrum_collection object from src to dest.
 *  must pass in a memory allocated SPECTRUM_COLLECTION_T* dest
 */
void copy_spectrum_collection(
  SPECTRUM_COLLECTION_T* src,///< spectrum to copy from -in
  SPECTRUM_COLLECTION_T* dest///< spectrum to copy to -out
  )
{
  SPECTRUM_T* new_spectrum;
  // copy each varible
  set_spectrum_collection_filename(dest,src->filename);
  set_spectrum_collection_comment(dest,src->comment);
  dest->num_charged_spectra = src->num_charged_spectra;
  dest->is_parsed = src->is_parsed;
  
  // copy spectrum
  SPECTRUM_ITERATOR_T* spectrum_iterator = new_spectrum_iterator(src);
  while(spectrum_iterator_has_next(spectrum_iterator)){
    new_spectrum = allocate_spectrum();
    copy_spectrum(spectrum_iterator_next(spectrum_iterator), new_spectrum);
    add_spectrum_to_end(dest, new_spectrum);
  }
  free_spectrum_iterator(spectrum_iterator);
}


// TESTME might have to check for bad format
/**
 * parses 'H' line into the spectrum_collection comments
 * all reminding comments are ignored if max length of comments are reached
 */
void parse_header_line(SPECTRUM_COLLECTION_T* spectrum_collection, FILE* file){
  long file_index = ftell(file); // stores the location of the current working line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0; // TODO maybe should be more sensible length?
  int new_line_length;
  int comment_field_length;

  while( (line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == 'H'){
      new_line_length = strlen(new_line);
      comment_field_length = strlen(spectrum_collection->comment);
      
      // check if max capacifty is too full for new comment
      if(new_line_length + comment_field_length + 1 < MAX_COMMENT){
        strcat(spectrum_collection->comment, new_line);
      }
      else{
        break; // exceed max comments
      }
    }
    else if(new_line[0] == 'S'){
      break; // end of 'H' lines
    }
    file_index = ftell(file);
  }
  free(new_line); // might change
  fseek(file, file_index, SEEK_SET);
}

/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
BOOLEAN_T parse_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< empty spectrum to parse into -out
)
{
  // spectrum_collection has already been parsed
  if(spectrum_collection->is_parsed){
    return FALSE;
  }

  FILE* file;
  SPECTRUM_T* parsed_spectrum;

  // check if file is still avaliable
  if ((file = fopen(spectrum_collection->filename,"r")) == NULL) {
    carp(CARP_ERROR, "File %s could not be opened",spectrum_collection->filename);
    return (FALSE);
  }
  // parse header lines 'H' into spectrum_collection comment 
  parse_header_line(spectrum_collection, file);

  //check to see if the mstoolkit is going to used.
  if (get_boolean_parameter("use-mstoolkit")) {
    //We now know that the file exists,
    //MSToolkit doesn't check or doesn't report an
    //error.  So we need this check about, but not
    //the file, so close it. (SJM)
    //also, I want to use parse_header_line for
    //the collection.
    fclose(file);
    carp(CARP_INFO,"using mstoolkit to parse spectra");
    MST_MSREADER_T* mst_reader = newMST_MSReader();
    MST_SPECTRUM_T* mst_spectrum = newMST_Spectrum();

    MST_MSReader_readFile(mst_reader, spectrum_collection->filename, mst_spectrum);

    while(MST_Spectrum_getScanNumber(mst_spectrum) != 0) {
      parsed_spectrum = allocate_spectrum();
      parse_spectrum_spectrum(parsed_spectrum, 
        mst_spectrum, 
        spectrum_collection->filename);
      if (!add_spectrum_to_end(spectrum_collection, parsed_spectrum)) {
        free_spectrum(parsed_spectrum);
        return FALSE;
      }
      MST_MSReader_readFile(mst_reader, NULL, mst_spectrum);
    }
    freeMST_Spectrum(mst_spectrum);
    freeMST_MSReader(mst_reader);
  } else {
    parsed_spectrum = allocate_spectrum();
    // parse one spectrum at a time
    while(parse_spectrum_file(parsed_spectrum, file, spectrum_collection->filename)){
      // is spectrum capacity not full?
      if(!add_spectrum_to_end(spectrum_collection, parsed_spectrum)){
        free_spectrum(parsed_spectrum);
        fclose(file);
        return FALSE;
      }
      parsed_spectrum = allocate_spectrum();
    }
    
    free_spectrum(parsed_spectrum); // CHECKME why free_spectrum??
    fclose(file);
  }
  spectrum_collection->is_parsed = TRUE;
  
  // good job!
  return TRUE;
}

// CHECKME test if capacity is correct might be off by one
/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum to the end of the spectra array
 * should only be used when the adding in increasing scan num order
 * when adding in random order should use add_spectrum
 * spectrum must be heap allocated
 */
BOOLEAN_T add_spectrum_to_end(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to add to spectrum_collection -in
  )
{
  // FIXME eventually might want it to grow dynamically
  // check if spectrum capacity is full
  if(get_spectrum_collection_num_spectra(spectrum_collection) == MAX_SPECTRA){
    carp(CARP_ERROR,"ERROR: cannot add spectrum, capacity full\n"); 
    return FALSE;
  }
  // set spectrum
  spectrum_collection->spectra[spectrum_collection->num_spectra] = spectrum;
  ++spectrum_collection->num_spectra;
  spectrum_collection->num_charged_spectra += get_spectrum_num_possible_z(spectrum);
  return TRUE;
}

/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum in correct order into the spectra array
 * spectrum must be heap allocated
 */
BOOLEAN_T add_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to add to spectrum_collection -in
  )
{
  // FIXME eventually might want it to grow dynamically
  int add_index = 0;
  int spectrum_index;
  
  // check if spectrum capacity is full
  if(get_spectrum_collection_num_spectra(spectrum_collection) == MAX_SPECTRA){
    carp(CARP_ERROR,"ERROR: cannot add spectrum, capacity full\n"); 
    return FALSE;
  }
  // find correct location
  for(; add_index < spectrum_collection->num_spectra; ++add_index){
    if(get_spectrum_first_scan(spectrum_collection->spectra[add_index])>
       get_spectrum_first_scan(spectrum)){
      break;
    }
  }
  // do we add to end?
  if(add_index != spectrum_collection->num_spectra +1){
    spectrum_index = spectrum_collection->num_spectra;
    // shift all spectrum that have greater or equal index to add_index to right  
    for(; spectrum_index >= add_index; --spectrum_index){
      spectrum_collection->spectra[spectrum_index+1] = 
        spectrum_collection->spectra[spectrum_index];
    }
  }
  
  // set spectrum
  spectrum_collection->spectra[add_index] = spectrum;
  ++spectrum_collection->num_spectra;
  spectrum_collection->num_charged_spectra += get_spectrum_num_possible_z(spectrum);
  return TRUE;
}


// FIXME maybe a faster way? can't perform binary search since we must know the array index
/**
 * Removes a spectrum from the spectrum_collection.
 */
void remove_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to be removed from spectrum_collection -in
  )
{
  int scan_num = get_spectrum_first_scan(spectrum);
  int spectrum_index = 0;
  
  // find where the spectrum is located in the spectrum array
  for(; spectrum_index < spectrum_collection->num_spectra; ++spectrum_index){
    if(scan_num ==
       get_spectrum_first_scan(spectrum_collection->spectra[spectrum_index])){
      break;
    }
  }
  
  free_spectrum(spectrum_collection->spectra[spectrum_index]);
  
  // shift all the spectra to the left to fill in the gap
  for(; spectrum_index < spectrum_collection->num_spectra; ++spectrum_index){
    spectrum_collection->spectra[spectrum_index] =
      spectrum_collection->spectra[spectrum_index+1];
  }
  
  --spectrum_collection->num_spectra;
  spectrum_collection->num_charged_spectra -= get_spectrum_num_possible_z(spectrum);
} 


/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan. Use binary search
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
BOOLEAN_T get_spectrum_collection_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< The spectrum collection -out
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  SPECTRUM_T* spectrum ///< The (empty) allocated SPECTRUM_T object -in
  )
{
  FILE* file;
  long target_index;
  // check if file is still avaliable
  if ((file = fopen(spectrum_collection->filename,"r")) == NULL) {
    carp(CARP_ERROR,"File %s could not be opened",spectrum_collection->filename);
    return (FALSE);
  }

  if (get_boolean_parameter("use-mstoolkit")) {
    //We now know that the file exists,
    //MSToolkit doesn't check or doesn't report an
    //error.  So we need this check about, but not
    //the file, so close it. (SJM)
    fclose(file);
    carp(CARP_INFO,"using mstoolkit to parse spectrum");
    MST_MSREADER_T* mst_reader = newMST_MSReader();
    MST_SPECTRUM_T* mst_spectrum = newMST_Spectrum();
    BOOLEAN_T parsed = FALSE;


    MST_MSReader_readFileScan(mst_reader, 
      spectrum_collection->filename, 
      mst_spectrum, first_scan);

    if (MST_Spectrum_getScanNumber(mst_spectrum) != 0) {
      parse_spectrum_spectrum(spectrum, 
            mst_spectrum, 
            spectrum_collection->filename);
      parsed = TRUE;
    }
    else {
      carp(CARP_ERROR,"Spectrum %d does not exist in file", first_scan);
      parsed = FALSE;
    }

    freeMST_Spectrum(mst_spectrum);
    freeMST_MSReader(mst_reader);
    return parsed;
  } else {

    target_index = binary_search_spectrum(file, first_scan);
    // first_scan not found
    if(target_index == -1){
      fclose(file);
      return FALSE;
    }
    fseek(file, target_index, SEEK_SET);
    // parse spectrum, check if failed to parse spectrum return false
    if(!parse_spectrum_file(spectrum, file, spectrum_collection->filename)){
      fclose(file);
      return FALSE;
    }
    fclose(file);
    return TRUE;
  }
}

/**
 * 
 *\returns the file position of the query spectrum begins in file
 * represented as an long integer, the return value of ftell()
 *\returns -1 if failed to find the query spectrum
 */
long binary_search_spectrum(
  FILE* file,///< the file to search -in
  int first_scan///< query scan num -in
)
{   
  long low_index = ftell(file); 
  long high_index;
  long mid_index;
  long working_index;
  long end_of_file_index;

  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  
  // set initial high and low position
  if(fseek(file, 1,SEEK_END) != -1){
    high_index = ftell(file);
    end_of_file_index = high_index;
  }
  else{
    carp(CARP_ERROR, "error: file corrupted");
    return -1;
  }  
  while(low_index <= high_index){

    mid_index = (low_index + high_index)/2;
    
    if(fseek(file, mid_index,SEEK_SET) != -1){
      
      working_index = ftell(file);
      // check each line until reach 'S' line
      while((line_length =  getline(&new_line, &buf_length, file)) != -1){

        if(new_line[0] == 'S'){
          int compare = 
            match_first_scan_line(new_line, buf_length, first_scan);
          if(compare == 0){
            free(new_line);
            return working_index; // found the query match
          }
          else if(compare == 1){
            high_index = mid_index-1;
            break;
          }
          else{ // compare == -1
            low_index = mid_index+1;
            break;
          }
        }
        // store the next working line
        working_index = ftell(file);
        
        // if the High is at EOF or working_index starts looking bellow high
        if(working_index > high_index || 
           working_index == end_of_file_index-1 ){ 
          high_index = mid_index-1;
          break;
        } 
      }
    }
    else{
      carp(CARP_ERROR, "error: file corrupted");
      return -1;
    }
  }
  free(new_line);
  return -1; // failed to find the query spectrum
}


/**
 * Parses the 'S' line of the a spectrum,
 * check if the first_scan match the query scan number
 * \returns 0 if match, 
 * \returns 1 if query_first_scan < first_scan
 * \returns -1 if query_first_scan > first_scan 
*/
int match_first_scan_line(
  char* line, ///< the line of interest -in
  int buf_length, ///< line length -in
  int query_first_scan ///< the query first scan num -in
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
    carp(
      CARP_FATAL,
      "Failed to parse 'S' line:\n %s "
      "Incorrect file format\n",
      line
    );
  }
  if(first_scan == query_first_scan){
    return 0;
  }
  else if(first_scan >  query_first_scan){
    return 1;
  }
  return -1; // first_scan <  query_first_scan
}

/******************************************************************************/

/**  ////// TESTME////
 * sets the filename of the ms2 file the spectra were parsed
 * this function should be used only the first time the filename is set
 * to change existing filename use set_spectrum_collection_filename
 * copies the value from arguement char* filename into a heap allocated memory
 */
void set_spectrum_collection_new_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save filename -out
  char* filename ///< filename -in
  )
{
  int filename_length = strlen(filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  spectrum_collection->filename =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * sets the filename of the ms2 file the spectrum_collection was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void set_spectrum_collection_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save filename -out
  char* filename ///< filename -in
  )
{
  free(spectrum_collection->filename);
  set_spectrum_collection_new_filename(spectrum_collection, filename);
}

/**
 * \returns the filename of the ms2 file the spectra was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* get_spectrum_collection_filename(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum collection's filename -in 
  )
{  
  int filename_length = strlen(spectrum_collection->filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename, spectrum_collection->filename, filename_length);  
}

/**
 * \returns the current number of spectrum in the spectrum_collection
 */
int get_spectrum_collection_num_spectra(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection save filename -in                                         
  )
{
  return spectrum_collection->num_spectra;
}

/**
 * \returns The current number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra) in the spectrum_collection
 */
int get_spectrum_collection_num_charged_spectra(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection save filename -in
  )
{
  return spectrum_collection->num_charged_spectra;
}


/**
 * \returns the comments from the spectrum_collection
 * the return char* points to a newly heap allocated copy of the comments
 * user must free the new string object
 */

char* get_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
  )
{
  char* comments = (char *)mycalloc(1, sizeof(char)*MAX_COMMENT);
  return strncpy(comments, spectrum_collection->comment, MAX_COMMENT); 
}

/**
 * sets the comment of the spectrum_collection
 * copies the new_comment into a newly heap allocated copy of the comment
 */

void set_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save comment -in                                         
  char* new_comment ///< the new comments to be copied
  )
{
  // is there enough memory for new comments?
  if(strlen(new_comment) + strlen(spectrum_collection->comment) +1 < MAX_COMMENT){
    strncat(spectrum_collection->comment, new_comment, MAX_COMMENT); 
  }
  else{
    carp(CARP_ERROR,"max comment exceeded\n");
  }
}


/**
 * \returns TRUE if the spectrum_collection file has been parsed
 */
BOOLEAN_T get_spectrum_collection_is_parsed(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
)
{
  return spectrum_collection->is_parsed;
}

/**
 * Takes the spectrum file name and creates a file with unique filenames.
 * The method will create one file for PSM result serializations for the 
 * target sequence and #number_decoy_set number of files for decoy PSM 
 * result serialization.  Thus, the FILE* array will contain,
 * at index 0, the target file and the allowed indices the decoy files.
 *
 * Template: "fileName_XXXXXX", where XXXXXX is random generated to be unique.
 * Also, sets psm_result_filenames pointer to the array of filenames for 
 * both the target and decoy psm results. The array is heap allocated, 
 * thus user must free it. Size is number_decoy_set +1 (for target)
 * \returns file handle array to the newly created files (target & decoy) 
 * and sets psm_result_filename.
 */
FILE** get_spectrum_collection_psm_result_filenames(
  SPECTRUM_COLLECTION_T* spectrum_collection, 
    ///< the spectrum_collection -in
  char* psm_result_folder_name, 
    ///< the folder name for where the result file should be placed -in
  char*** psm_result_filenames, 
    ///< pointer to be set to the array of filenames -out
  int number_decoy_set,  
    ///< the number of decoy sets to produce -in
  char* file_extension 
    ///< the file extension of the spectrum file (i.e. ".ms2") -in
  )
{
  int file_descriptor = -1;
  char suffix[25];
  // total number of files to create, target plus how many decoys needed
  int total_files = number_decoy_set + 1; 

  // check if psm_result_folder exist?
  if(access(psm_result_folder_name, F_OK)){
    // create PSM result folder
    if(mkdir(psm_result_folder_name, S_IRWXU+S_IRWXG+S_IRWXO) != 0){
      carp(CARP_ERROR, "failed to create psm result folder: %s", 
          psm_result_folder_name);
    }
  }
  
  // create FILE* array
  FILE** file_handle_array = (FILE**)mycalloc(total_files, sizeof(FILE*));
  char** filename_array = (char**)mycalloc(total_files, sizeof(char*));
  
  // extract filename from absolute path
  char** spectrum_file_path = 
    parse_filename_path(spectrum_collection->filename);
  
  // create a filename template that has psm_result_folder_name/spectrum_fname
  char* filename_template 
    = get_full_filename(psm_result_folder_name, spectrum_file_path[0]); 
  
  // now create files for first target file and then for decoys
  int file_idx;
  for(file_idx = 0; file_idx < total_files; ++file_idx){

    // is it target?
    if(file_idx == 0){

      // generate psm_result filename as 
      // psm_result_folder_name/spectrum_filename_XXXXXX
      filename_array[file_idx] = generate_name(filename_template, "_XXXXXX", 
          file_extension, "crux_match_target_");
    }

    // for decoys
    else{
      sprintf(suffix, "crux_match_decoy_%d_", file_idx);
      filename_array[file_idx] = generate_name(filename_template, "_XXXXXX", 
          file_extension, suffix);
    }

    // now open file handle
    if((file_descriptor = mkstemp(filename_array[file_idx])) == -1 ||
       (file_handle_array[file_idx] = fdopen(file_descriptor, "w+")) == NULL){
      
      // did we successfully create a file?
      if(file_descriptor != -1){
        unlink(filename_array[file_idx]);
        close(file_descriptor);
      }
      free(spectrum_file_path[0]);
      free(spectrum_file_path[1]);
      free(spectrum_file_path);
      free(filename_template);
      carp(CARP_ERROR, "failed to create PSM output file");
      return NULL;
    }
    // set permission for the file
    chmod(filename_array[file_idx], 0664);
  }
  
  free(spectrum_file_path[0]);
  free(spectrum_file_path[1]);
  free(spectrum_file_path);
  free(filename_template);
  
  // set output for result filenames
  *psm_result_filenames = filename_array;
  return file_handle_array;
}

/**
 * <int: number spectra> <--this will be over written by serialize_total_number_of_spectra method
 * <int: number of spectrum features>
 * <int: number of top ranked peptides serialized per spectra>
 *
 * 
 * // FIXME <int: ms2 file length><char*: ms2 filename>
 * // FIXME <int: fasta file length><char*: fasta filename>
 *
 * Serializes the header information for the binary PSM serialized files
 * Must run in pair with serialize_total_number_of_spectra.
 *
 * General order is, 
 * serialize_header -> serialize_psm_features -> serialize_total_number_of_spectra
 *\returns TRUE if serialized header successfully, else FALSE
 */
BOOLEAN_T serialize_header(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection -in
  char* fasta_file, ///< the fasta file 
  FILE* psm_file ///< the file to serialize the header information -out
  )
{
  int num_spectrum_features = 0;
  // set max number of matches to be serialized per spectrum
  int number_top_rank_peptide = get_int_parameter("top-match");
  char* file_fasta = parse_filename(fasta_file);
  // int file_fasta_length = strlen(file_fasta);
  char* file_ms2 = parse_filename(spectrum_collection->filename);
  // int file_ms2_length = strlen(file_ms2);
  
  // FIXME later if you want to be selective on charge to run
  // must refine the current serialize methods
  
  // num_charged_spectra will be over written by serialize_total_number_of_spectra method
  fwrite(&(spectrum_collection->num_charged_spectra), sizeof(int), 1, psm_file);
  fwrite(&(num_spectrum_features), sizeof(int), 1, psm_file);
  fwrite(&(number_top_rank_peptide), sizeof(int), 1, psm_file);
  
  carp(CARP_DETAILED_DEBUG, "Serialize header wrote %i top matches", 
       number_top_rank_peptide);
  // free up files
  free(file_ms2);
  free(file_fasta);
  
  return TRUE;
}

/**
 * Modifies the serialized header information for the binary PSM serialized files
 * Sets the total number of spectra serialized in the file
 * Assumes the first field in the file is the number total spectra serialized
 * Must be run after serialize_header
 *
 * General order is, 
 * serialize_header -> serialize_psm_features -> serialize_total_number_of_spectra
 *
 *\returns TRUE if total number of spectra seerialized in the file, else FALSE
 */
BOOLEAN_T serialize_total_number_of_spectra(
  int spectra_idx, ///< the number of spectra serialized in PSM file -in 
  FILE* psm_file ///< the file to serialize the header information -out
  )
{
  if( psm_file == NULL ){
    return FALSE;
  }
  // set to begining of file
  rewind(psm_file);

  // serialize the total number of spectra seerialized in the file
  fwrite(&(spectra_idx), sizeof(int), 1, psm_file);

  return TRUE;
}
/******************************************************************************/

/**
 * Instantiates a new spectrum_iterator object from spectrum_collection.
 * \returns a SPECTRUM_ITERATOR_T object.
 */
SPECTRUM_ITERATOR_T* new_spectrum_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectrum_collection to iterate -in
  )
{
  SPECTRUM_ITERATOR_T* spectrum_iterator =
    (SPECTRUM_ITERATOR_T*)mycalloc(1, sizeof(SPECTRUM_ITERATOR_T));
  spectrum_iterator->spectrum_collection = spectrum_collection;
  spectrum_iterator->spectrum_index = 0;
  return spectrum_iterator;
}
        
/**
 * Instantiates a new spectrum_iterator object from
 * spectrum_collection.  This iterator returns unique spectrum-charge
 * pairs (e.g.a spectrum to be searched as +2 and +3 is returned once
 * as +2 and once as +3).  The charge is returned by setting the int
 * pointer in the argument list.  The iterator also filters spectra by
 * mass so that none outside the spectrum-min-mass--spectrum-max-mass
 * range (as defined in parameter.c).
 * \returns a SPECTRUM_ITERATOR_T object.
 */
FILTERED_SPECTRUM_CHARGE_ITERATOR_T* new_filtered_spectrum_charge_iterator(
  SPECTRUM_COLLECTION_T* spectrum_collection///< spectra to iterate over
){        
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator =
  (FILTERED_SPECTRUM_CHARGE_ITERATOR_T*)mycalloc(1, 
                               sizeof(FILTERED_SPECTRUM_CHARGE_ITERATOR_T));
  iterator->spectrum_collection = spectrum_collection;
  iterator->has_next = FALSE;
  iterator->spectrum_index = -1;
  iterator->charges = NULL;
  iterator->num_charges = 0;
  iterator->charge_index = -1;
  iterator->min_mz = get_double_parameter("spectrum-min-mass");
  iterator->max_mz = get_double_parameter("spectrum-max-mass");
  const char* charge_str = get_string_parameter_pointer("spectrum-charge");
  if( strcmp( charge_str, "all") == 0){
    iterator->search_charge = 0;
  }else{
    iterator->search_charge = atoi(charge_str);
  }
  // queue next spectrum
  queue_next_spectrum(iterator);
  return iterator;
}

/**
 * Frees an allocated spectrum_iterator object.
 */
void free_spectrum_iterator(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< free spectrum_iterator -in
)
{
  free(spectrum_iterator);
}

/**
 * Frees an filtered_spectrum_charge_iterator object.
 */
void free_filtered_spectrum_charge_iterator(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator///< free spectrum_iterator -in
)
{
  free(iterator);
}

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T spectrum_iterator_has_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< is there a next spectrum? -in
)
{
  return (spectrum_iterator->spectrum_index 
          < get_spectrum_collection_num_spectra(
                    spectrum_iterator->spectrum_collection));
}

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T filtered_spectrum_charge_iterator_has_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator){

  if( iterator == NULL ){
    carp(CARP_ERROR, "Cannot find next for NULL iterator");
    return FALSE;
  }
  return iterator->has_next;
}

/**
 * The basic iterator function next.
 */
SPECTRUM_T* spectrum_iterator_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< return the next spectrum -in
)
{
  SPECTRUM_T* next_spectrum = 
    spectrum_iterator->spectrum_collection->spectra[spectrum_iterator->spectrum_index];
  ++spectrum_iterator->spectrum_index;
  return next_spectrum;
}

/**
 * The basic iterator function next.  Also returns the charge state to
 * use for this spectrum.
 */
SPECTRUM_T* filtered_spectrum_charge_iterator_next(
  FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator,
  int* charge                            ///< put charge here -out
)
{
  SPECTRUM_T* next_spectrum = 
    iterator->spectrum_collection->spectra[iterator->spectrum_index];
  *charge = iterator->charges[iterator->charge_index];

  queue_next_spectrum(iterator);

  return next_spectrum;
}

/**
 * \brief Sets up an iterator with the next spectrum that complies
 * with the constraints.  Sets has_next to FALSE when there are no
 * more spectra in the collection that pass.  Increments
 * spectrum_index and charge_index.
 */
void queue_next_spectrum(FILTERED_SPECTRUM_CHARGE_ITERATOR_T* iterator){
  if( iterator == NULL ){
    carp(CARP_ERROR, "Cannot queue spectrum for NULL iterator.");
    return;
  }

  SPECTRUM_T* spec = NULL;

  // Are there any more charge states for this spectrum?
  if( iterator->charge_index < iterator->num_charges-1 ){
    iterator->charge_index++;
  }
  // Are there any more spectra?
  else if( iterator->spectrum_index < get_spectrum_collection_num_spectra(
                                 iterator->spectrum_collection)-1 ){
    iterator->spectrum_index++;
    spec = iterator->spectrum_collection->spectra[iterator->spectrum_index];
    // first free any existing charges in the iterator
    if( iterator->charges ){
      free(iterator->charges);
    }
    iterator->num_charges = get_charges_to_search(spec, &(iterator->charges));
    iterator->charge_index = 0;
  }else{ // none left
    iterator->has_next = FALSE;
    return;
  }

  // Does the current pass?
  spec = iterator->spectrum_collection->spectra[iterator->spectrum_index];
  int this_charge = -1;
  if (iterator->charge_index < iterator->num_charges) {
    this_charge = iterator->charges[iterator->charge_index];
  }
  double mz = get_spectrum_precursor_mz(spec);

  if( iterator->search_charge == 0 || iterator->search_charge == this_charge ){
    if( mz >= iterator->min_mz && mz <= iterator->max_mz ){
      // passes both tests
      iterator->has_next = TRUE;
      return;
    }
  }
  
  // try the next spectrum
  queue_next_spectrum(iterator);

}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
