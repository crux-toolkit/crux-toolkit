/**
 * \file MS2SpectrumCollection.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * \brief Class to read spectra from .ms2 files.
 */
#include "MS2SpectrumCollection.h" 
#include "crux-utils.h"
#include "parameter.h"
#include "Spectrum.h"

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does not parse file. 
 */
MS2SpectrumCollection::MS2SpectrumCollection(
  const char* filename   ///< The spectrum collection filename.
 ) : SpectrumCollection(filename){
   memset(comment_, 0, MAX_COMMENT);
}

/**
 * Parses all 'H' line into the spectrum_collection comments until the
 * maximum size of the comment buffer has been reached.
 */
void MS2SpectrumCollection::parseHeaderLine(
  FILE* file
){

  long file_index = ftell(file); // location of the current line in the file
  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0; // TODO maybe should be more sensible length?
  unsigned int new_line_length;
  unsigned int comment_field_length;

  while( (line_length =  getline(&new_line, &buf_length, file)) != -1){
    if(new_line[0] == 'H'){
      new_line_length = strlen(new_line);
      comment_field_length = strlen(comment_);
      
      // check if max capacifty is too full for new comment
      if(new_line_length + comment_field_length + 1 < MAX_COMMENT){
        strcat(comment_, new_line);
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
 * \returns True if the spectra are parsed successfully. False if otherwise.
 */
bool MS2SpectrumCollection::parse() {

  // spectrum_collection has already been parsed
  if(is_parsed_){
    return false;
  }

  FILE* file;

  // get a list of scans to include if requested
  const char* range_string = get_string_parameter("scan-number");
  int first_scan;
  int last_scan;
  bool success = get_range_from_string(range_string, first_scan, last_scan);

  if( !success ){
    carp(CARP_FATAL, "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.", range_string);
  }
  
  // check if file is still avaliable
  if ((file = fopen(filename_.c_str(),"rb")) == NULL) {
    carp(CARP_ERROR, "Spectrum file %s could not be opened",
         filename_.c_str());
    return false;
  }

  parseHeaderLine(file);

  // parse one spectrum at a time
  Spectrum* parsed_spectrum = Spectrum::newSpectrumMs2(file, filename_.c_str());
  while(parsed_spectrum){
    // include this a scan?
    if( parsed_spectrum->getFirstScan() < first_scan ){
      delete parsed_spectrum;
      parsed_spectrum = Spectrum::newSpectrumMs2(file, filename_.c_str());
      continue;
    } 
    // are we past the last scan?
    if( parsed_spectrum->getFirstScan() > last_scan ){
      break;
    }
    this->addSpectrumToEnd(parsed_spectrum);
    parsed_spectrum = Spectrum::newSpectrumMs2(file, filename_.c_str());
  }
  //    delete parsed_spectrum; // CHECKME why free_spectrum??
  fclose(file);
  is_parsed_ = true;
  return true;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.  Use binary search.  Removes any
 * existing information in the given spectrum.
 * \returns True if the spectrum was allocated, false on error.
 */
bool MS2SpectrumCollection::getSpectrum(
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  Spectrum* spectrum   ///< Put the spectrum info here
  )
{
  FILE* file;
  long target_index;
  // check if file is still avaliable
  if ((file = fopen(filename_.c_str(), "rb")) == NULL) {
    carp(CARP_ERROR, "File %s could not be opened", filename_.c_str());
    return false;
  }

  if( spectrum == NULL ){
    carp(CARP_ERROR, "Can't parse into a NULL spectrum.");
    return false;
  }

  target_index = binarySearchSpectrum(file, first_scan);
  // first_scan not found
  if(target_index == -1){
    fclose(file);
    return false;
  }
  fseek(file, target_index, SEEK_SET);
  // parse spectrum, check if failed to parse spectrum return false
  if(!spectrum->parseMs2(file, filename_.c_str())){
    fclose(file);
    return false;
  }
  fclose(file);
  return true;
}

/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan. Use binary search
 * \returns The spectrum data from file or NULL.
 */
Spectrum* MS2SpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
  )
{
  FILE* file;
  long target_index;
  // check if file is still avaliable
  if ((file = fopen(filename_.c_str(), "rb")) == NULL) {
    carp(CARP_ERROR, "File %s could not be opened", filename_.c_str());
    return (false);
  }

  Spectrum* return_spec = NULL;
  target_index = binarySearchSpectrum(file, first_scan);
  // first_scan not found
  if(target_index == -1){
    fclose(file);
    return NULL;
  }
  fseek(file, target_index, SEEK_SET);
  // parse spectrum, check if failed to parse spectrum return false
  return_spec = Spectrum::newSpectrumMs2(file, filename_.c_str());
  fclose(file);
  return return_spec;
}

/**
 * \returns The file position of the query spectrum in the file,
 * i.e. the return value of ftell().  Or -1 if failed to find the
 * spectrum.
 */
long MS2SpectrumCollection::binarySearchSpectrum(
  FILE* file,    ///< the file to search 
  int first_scan ///< query scan num 
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
            matchFirstScanLine(new_line, buf_length, first_scan);
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
int MS2SpectrumCollection::matchFirstScanLine(
  char* line,          ///< the line of interest
  int buf_length,      ///< line length 
  int query_first_scan ///< the query first scan num
  )
{
  
  char* spliced_line = new char[buf_length];
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

  delete [] spliced_line ;

  if(first_scan == query_first_scan){
    return 0;
  }
  else if(first_scan >  query_first_scan){
    return 1;
  }
  return -1; // first_scan <  query_first_scan
}

/**
 * \returns The comments from the .ms2 file.
 */
const char* MS2SpectrumCollection::getComment()
{
  return comment_;
}

/**
 * Sets or adds to the comment of the ms2 file.
 */

void MS2SpectrumCollection::setComment(
  const char* new_comment ///< the new comments to be copied
  )
{
  // is there enough memory for new comments?
  if(strlen(new_comment) + strlen(comment_) +1 < MAX_COMMENT){
    strncat(comment_, new_comment, MAX_COMMENT); 
  }
  else{
    carp(CARP_ERROR,"max comment exceeded\n");
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
