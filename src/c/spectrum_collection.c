/*****************************************************************************
 * \file spectrum_collection.c
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * DESCRIPTION: code to support working with collection of multiple spectra
 * REVISION: $Revision: 1.12 $
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

#define MAX_SPECTRA 40000 ///< max number of spectrums



long binary_search_spectrum(FILE* file, int first_scan);

int match_first_scan_line(
  char* line, 
  int buf_length, 
  int query_first_scan
  );

/**
 * \struct spectrum_collection 
 * Holds one or more spectrum objects.
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
 * Allows iteration over the spectra within a spectrum_collection.
 */
struct spectrum_iterator {
  SPECTRUM_COLLECTION_T* spectrum_collection;///< The collection whose spectrum to iterate over. 
  int  spectrum_index; ///< The index of the current spectrum;
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
  char* filename///< The spectrum collection filename. -in
  )
{
  SPECTRUM_COLLECTION_T* spectrum_collection =  allocate_spectrum_collection();

  if(access(filename, F_OK)){
    fprintf(stderr,"File %s could not be opened\n",filename);
    free(spectrum_collection);
    exit(1);
  } //FIXME check if file is empty
  
  set_spectrum_collection_filename(spectrum_collection, filename); ///FIXME might need to use new
  
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
  free(spectrum_collection->comment);
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
  
  fprintf(file,"comment: %s\n", spectrum_collection->comment);
  fprintf(file,"filename: %s\n", spectrum_collection->filename);
  
  //print each spectrum
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
  //copy each varible
  set_spectrum_collection_filename(dest,src->filename);
  //set_spectrum_collection_comment(dest,src->comment);
  dest->is_parsed = src->is_parsed;
  
  //copy spectrum
  SPECTRUM_ITERATOR_T* spectrum_iterator = new_spectrum_iterator(src);
  while(spectrum_iterator_has_next(spectrum_iterator)){
    add_spectrum(dest, spectrum_iterator_next(spectrum_iterator));
  }
  free_spectrum_iterator(spectrum_iterator);
}

//FIXME must be able to parse 'H' line
/**
 * Parses all the spectra from file designated by the filename member
 * variable.
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
BOOLEAN_T parse_spectrum_collection(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< empty spectrum to parse into -out
)
{
  //spectrum_collection has already been parsed
  if(spectrum_collection->is_parsed){
    return FALSE;
  }

  FILE* file;
  SPECTRUM_T* parsed_spectrum;

  //check if file is still avaliable
  if ((file = fopen(spectrum_collection->filename,"r")) == NULL) {
    fprintf(stderr,"File %s could not be opened",spectrum_collection->filename);
    return (FALSE);
  }
  
  parsed_spectrum = allocate_spectrum();
  //parse one spectrum at a time
  while(parse_spectrum_file(parsed_spectrum, file)){
    //is spectrum capacity not full?
    if(!add_spectrum(spectrum_collection, parsed_spectrum)){
      free_spectrum(parsed_spectrum);
      fclose(file);
      return FALSE;
    }
    parsed_spectrum = allocate_spectrum();
  }
  
  free_spectrum(parsed_spectrum); //CHECKME why free_spectrum??
  fclose(file);

  spectrum_collection->is_parsed = TRUE;
  
  //good job!
  return TRUE;
}

//CHECKME test if capacity is correct might be off by one
/**
 * Adds a spectrum to the spectrum_collection.
 * spectrum must be heap allocated
 */
BOOLEAN_T add_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to add to spectrum_collection -in
  )
{
  //FIXME eventually might want it to grow dynamically
  //check if spectrum capacity is full
  if(get_spectrum_collection_num_spectra(spectrum_collection) == MAX_SPECTRA){
    fprintf(stderr,"ERROR: cannot add spectrum, capacity full\n"); 
    return FALSE;
  }
  //set spectrum
  spectrum_collection->spectra[spectrum_collection->num_spectra] = spectrum;
  ++spectrum_collection->num_spectra;
  return TRUE;
}


/**
 * Removes a spectrum from the spectrum_collection.
 */
/*
void remove_spectrum(
  SPECTRUM_COLLECTION_T* spectrum_collection,///< the working spectrum_collection -out
  SPECTRUM_T* spectrum ///< spectrum to be removed from spectrum_collection -in
  )
{
  
  ///WRITEME

} 
*/

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
long binary_search_spectrum(
  FILE* file,///< the file to search -in
  int first_scan///< query scan num -in
)
{   
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

/******************************************************************************/

/**  //////TESTME////
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
  int filename_length = strlen(filename) +1; //+\0
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
  int filename_length = strlen(spectrum_collection->filename) +1; //+\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename,spectrum_collection->filename,filename_length);  
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
 * \returns the comments from the spectrum_collection
 * the return char* points to a newly heap allocated copy of the comments
 */
/*
char* get_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
  )
{
  return "yeah!";

}
*/

/**
 * sets the comment of the spectrum_collection
 * copies the new_comment into a newly heap allocated copy of the comment
 */
/*
void set_spectrum_collection_comment(
  SPECTRUM_COLLECTION_T* spectrum_collection, ///< the spectrum_collection save comment -in                                         
  char* new_comment ///< the new comments to be copied
  )
{
  ///WRITEME
}
*/

/**
 * \returns TRUE if the spectrum_collection file has been parsed
 */
BOOLEAN_T get_spectrum_collection_is_parsed(
  SPECTRUM_COLLECTION_T* spectrum_collection ///< the spectrum_collection -in                                         
)
{
  return spectrum_collection->is_parsed;
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
 * Frees an allocated spectrum_iterator object.
 */
void free_spectrum_iterator(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< free spectrum_iterator -in
)
{
  free(spectrum_iterator);
}

/**
 * The basic iterator function has_next.
 */
BOOLEAN_T spectrum_iterator_has_next(
  SPECTRUM_ITERATOR_T* spectrum_iterator///< is there a next spectrum? -in
)
{
  return (spectrum_iterator->spectrum_index 
          < get_spectrum_collection_num_spectra(spectrum_iterator->spectrum_collection));
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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
