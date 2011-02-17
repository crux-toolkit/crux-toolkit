/*************************************************************************//**
 * \file SpectrumCollection.cpp
 * AUTHOR: Chris Park
 * CREATE DATE: 28 June 2006
 * \brief code to support working with collection of multiple spectra
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <vector>
#include "objects.h"
#include "Spectrum.h"
#include "SpectrumCollection.h" 
#include "protein_index.h" 
#include "peak.h"
#include "utils.h"
#include "unistd.h"
#include "parameter.h"

#include "Spectrum.h"
#include "MSReader.h"

/**
 * \returns An (empty) heap allocated spectrum_collection object.
 */
SpectrumCollection::SpectrumCollection() {

  filename_ = NULL;
  memset(comment_, 0, sizeof(char)*MAX_COMMENT);
  is_parsed_ = false;
}

/**
 * Instantiates a new spectrum_collection object from a filename. 
 * Does *NOT* parse all spectra in the file. 
 * This will be done lazily depending on the subsequent method
 * calls (parse() or  getSpectrum()).
 * \returns  SPECTRUM_COLLECTION_T
 */
SpectrumCollection::SpectrumCollection(
  const char* filename///< The spectrum collection filename. -in
  ) {

  filename_ = NULL;
  memset(comment_, 0, sizeof(char)*MAX_COMMENT);
  is_parsed_ = false;
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
    free(absolute_path_file);
    carp(CARP_FATAL,"File %s could not be opened\n", absolute_path_file);
  } // FIXME check if file is empty
  this->setFilename(absolute_path_file);
  #ifndef DARWIN
  free(absolute_path_file);
  #endif
}

/**
 * Copy constructor.
 */
SpectrumCollection::SpectrumCollection(
  SpectrumCollection& old_collection
  ) {

  // copy each varible
  this->setFilename(old_collection.filename_);
  this->setComment(old_collection.comment_);
  is_parsed_ = old_collection.is_parsed_;
  num_charged_spectra_ = old_collection.num_charged_spectra_;
  // copy spectrum
  for (SpectrumIterator spectrum_iterator = old_collection.begin();
    spectrum_iterator != old_collection.end();
    ++spectrum_iterator) {

    Spectrum* old_spectrum = *spectrum_iterator;
    Spectrum* new_spectrum = new Spectrum(*old_spectrum);
    this->addSpectrumToEnd(new_spectrum);
  }
} 

/**
 * Default destructor
 */
SpectrumCollection::~SpectrumCollection() {

  for (SpectrumIterator spectrum_iterator = this->begin();
    spectrum_iterator != this->end();
    ++spectrum_iterator) {
    delete *spectrum_iterator;    
  }
  spectra_.clear();
  free(filename_);
}  

/**
 * \returns the begining of the spectra vector
 */
SpectrumIterator SpectrumCollection::begin() {
  return spectra_.begin();
}

/**
 * \returns the end of the spectra vector
 */
SpectrumIterator SpectrumCollection::end() {
  return spectra_.end();
}

/**
 * Prints a spectrum_collection object to file.
 */
void SpectrumCollection::printSpectrumCollection(
  FILE* file ///< file for output -out
  )
{
  
  fprintf(file,"comment:\n%s", comment_);
  fprintf(file,".ms2 Filename: %s\n", filename_);
  
  // print each spectrum
  for (SpectrumIterator spectrum_iterator = this->begin();
    spectrum_iterator != this->end();
    ++spectrum_iterator) {

    (*spectrum_iterator)->print(file);
  }
}

// TESTME might have to check for bad format
/**
 * parses 'H' line into the spectrum_collection comments
 * all reminding comments are ignored if max length of comments are reached
 */
void SpectrumCollection::parseHeaderLine(
  FILE* file
){

  long file_index = ftell(file); // stores the location of the current working line in the file
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
 * \returns TRUE if the spectra are parsed successfully. FALSE if otherwise.
 */
bool SpectrumCollection::parse() {

  // spectrum_collection has already been parsed
  if(is_parsed_){
    return false;
  }

  FILE* file;
  Spectrum* parsed_spectrum;

  // get a list of scans to include if requested
  const char* range_string = get_string_parameter("scan-number");
  int first_scan = get_first_in_range_string(range_string);
  int last_scan = get_last_in_range_string(range_string);
  if( first_scan == -1 || last_scan == -1 ){
    carp(CARP_FATAL, "The scan number range '%s' is invalid. "
         "Must be of the form <first>-<last>.", range_string);
  }
  
  // check if file is still avaliable
  if ((file = fopen(filename_,"r")) == NULL) {
    carp(CARP_ERROR, "File %s could not be opened",
         filename_);
    return false;
  }

  if (!get_boolean_parameter("use-mgf")) {
    // parse header lines 'H' into spectrum_collection comment 
    parseHeaderLine(file);

  }

  //check to see if the mstoolkit is going to used.
  if (get_boolean_parameter("use-mstoolkit")) {
    //We now know that the file exists,
    //MSToolkit doesn't check or doesn't report an
    //error.  So we need this check about, but not
    //the file, so close it. (SJM)
    //also, I want to use parse_header_line for
    //the collection.
    fclose(file);
    carp(CARP_INFO, "Using mstoolkit to parse spectra.");

    MSToolkit::MSReader* mst_reader = new MSToolkit::MSReader();
    MSToolkit::Spectrum* mst_spectrum = new MSToolkit::Spectrum();

    // only read ms2 scans
    mst_reader->setFilter(MSToolkit::MS2);


    mst_reader -> readFile(filename_, *mst_spectrum);

    while(mst_spectrum -> getScanNumber() != 0) {
      // is this a scan to include? if not skip it
      if( mst_spectrum -> getScanNumber() < first_scan ){
        mst_reader -> readFile(NULL, *mst_spectrum);
        continue;
      } 
      // are we past the last scan?
      if( mst_spectrum ->  getScanNumber() > last_scan ){
        break;
      }
      parsed_spectrum = new Spectrum();
      parsed_spectrum->parseMstoolkitSpectrum(mst_spectrum, filename_);

      if (!this->addSpectrumToEnd(parsed_spectrum)) {
        delete parsed_spectrum;
        return FALSE;
      }
      mst_reader -> readFile(NULL, *mst_spectrum);
    }
    delete mst_spectrum;
    delete mst_reader;
  } else { // not MSToolkit
    // parse one spectrum at a time
    Spectrum* parsed_spectrum = 
      Spectrum::newSpectrumFromFile(file, filename_);
    while(parsed_spectrum){
      // is this a scan to include? if not skip it
      if( parsed_spectrum->getFirstScan() < first_scan ){
        delete parsed_spectrum;
        parsed_spectrum = 
          Spectrum::newSpectrumFromFile(file, filename_);
        continue;
      } 
      // are we past the last scan?
      if( parsed_spectrum->getFirstScan() > last_scan ){
        break;
      }
      // is spectrum capacity not full?
      if(!this->addSpectrumToEnd(parsed_spectrum)){
        delete parsed_spectrum;
        fclose(file);
        return false;
      }
      parsed_spectrum = 
        Spectrum::newSpectrumFromFile(file, filename_);
    }
    
    delete parsed_spectrum; // CHECKME why free_spectrum??
    fclose(file);
  }
  is_parsed_ = TRUE;
  
  // good job!
  return true;
}

// CHECKME test if capacity is correct might be off by one
/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum to the end of the spectra array
 * should only be used when the adding in increasing scan num order
 * when adding in random order should use add_spectrum
 * spectrum must be heap allocated
 */
bool SpectrumCollection::addSpectrumToEnd(
  Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  )
{
  // set spectrum
  spectra_.push_back(spectrum);
  num_charged_spectra_ += spectrum->getNumZStates();
  return true;
}

/**
 * Adds a spectrum to the spectrum_collection.
 * adds the spectrum in correct order into the spectra array
 * spectrum must be heap allocated
 */
bool SpectrumCollection::addSpectrum(
  Spectrum* spectrum ///< spectrum to add to spectrum_collection -in
  )
{
  unsigned int add_index = 0;

  // find correct location
  // TODO -- replace with binary search if necessary.
  for(; add_index < spectra_.size(); ++add_index){
    if(spectra_[add_index]->getFirstScan() >
       spectrum->getFirstScan()){
      break;
    }
  }

  spectra_.insert(spectra_.begin()+add_index, spectrum);

  num_charged_spectra_ += spectrum->getNumZStates();
  return true;
}


// FIXME maybe a faster way? can't perform binary search since we must know the array index
/**
 * Removes a spectrum from the spectrum_collection.
 */
void SpectrumCollection::removeSpectrum(
  Spectrum* spectrum ///< spectrum to be removed from spectrum_collection -in
  )
{
  int scan_num = spectrum->getFirstScan();
  unsigned int spectrum_index = 0;
  
  // find where the spectrum is located in the spectrum array
  for(; spectrum_index < spectra_.size(); ++spectrum_index){
    if(scan_num == spectra_[spectrum_index]->getFirstScan() ){
      break;
    }
  }
  
  num_charged_spectra_ -= spectrum->getNumZStates();

  delete spectra_[spectrum_index];
  spectra_[spectrum_index] = NULL;
  spectra_.erase(spectra_.begin() + spectrum_index);

} 


/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan.  Use binary search.  Removes any
 * existing information in the given spectrum.
 * \returns TRUE if the spectrum was allocated, FALSE on error.
 */
bool SpectrumCollection::getSpectrum(
  int first_scan,      ///< The first scan of the spectrum to retrieve -in
  Spectrum* spectrum   ///< Put the spectrum info here
  )
{
  FILE* file;
  long target_index;
  // check if file is still avaliable
  if ((file = fopen(filename_,"r")) == NULL) {
    carp(CARP_ERROR,"File %s could not be opened",filename_);
    return (false);
  }

  if( spectrum == NULL ){
    carp(CARP_ERROR, "Can't parse into a NULL spectrum.");
    return false;
  }

  if (get_boolean_parameter("use-mstoolkit")) {
    //We now know that the file exists,
    //MSToolkit doesn't check or doesn't report an
    //error.  So we need this check about, but not
    //the file, so close it. (SJM)
    fclose(file);
    carp(CARP_INFO,"using mstoolkit to parse spectrum");
    MSToolkit::MSReader* mst_reader = new MSToolkit::MSReader();
    MSToolkit::Spectrum* mst_spectrum = new MSToolkit::Spectrum();
    bool parsed = false;

    mst_reader -> readFile(
      filename_,
      *mst_spectrum,
      first_scan);

    if (mst_spectrum -> getScanNumber() != 0) {
      spectrum->parseMstoolkitSpectrum(mst_spectrum,
                                         filename_);
      parsed = true;
    }
    else {
      carp(CARP_ERROR,"Spectrum %d does not exist in file", first_scan);
      parsed = false;
    }
    delete mst_spectrum;
    delete mst_reader;
    return parsed;
  } else {

    target_index = binarySearchSpectrum(file, first_scan);
    // first_scan not found
    if(target_index == -1){
      fclose(file);
      return false;
    }
    fseek(file, target_index, SEEK_SET);
    // parse spectrum, check if failed to parse spectrum return false
    if(!spectrum->parseFile(file, filename_)){
      fclose(file);
      return false;
    }
    fclose(file);
    return true;
  }
}


/**
 * Parses a single spectrum from a spectrum_collection with first scan
 * number equal to first_scan. Use binary search
 * \returns TRUE if the spectrum with. FALSE is failure.
 */
Spectrum* SpectrumCollection::getSpectrum(
  int first_scan      ///< The first scan of the spectrum to retrieve -in
  )
{
  FILE* file;
  long target_index;
  // check if file is still avaliable
  if ((file = fopen(filename_,"r")) == NULL) {
    carp(CARP_ERROR,"File %s could not be opened",filename_);
    return (false);
  }

  Spectrum* return_spec = NULL;

  if (get_boolean_parameter("use-mstoolkit")) {
    //We now know that the file exists,
    //MSToolkit doesn't check or doesn't report an
    //error.  So we need this check about, but not
    //the file, so close it. (SJM)
    fclose(file);
    carp(CARP_INFO,"using mstoolkit to parse spectrum");
    MSToolkit::MSReader* mst_reader = new MSToolkit::MSReader();
    MSToolkit::Spectrum* mst_spectrum = new MSToolkit::Spectrum();

    mst_reader -> readFile(
      filename_,
      *mst_spectrum,
      first_scan);

    if (mst_spectrum -> getScanNumber() != 0) {
      return_spec = new Spectrum();
      return_spec->parseMstoolkitSpectrum(mst_spectrum,
                                            filename_); 
    } else {
      carp(CARP_ERROR,"Spectrum %d does not exist in file", first_scan);
      return_spec = NULL;
    }
    delete mst_spectrum;
    delete mst_reader;
  } else {

    target_index = binarySearchSpectrum(file, first_scan);
    // first_scan not found
    if(target_index == -1){
      fclose(file);
      return NULL;
    }
    fseek(file, target_index, SEEK_SET);
    // parse spectrum, check if failed to parse spectrum return false
    return_spec = 
      Spectrum::newSpectrumFromFile(file, filename_);
    fclose(file);
  }
  return return_spec;
}

/**
 * 
 *\returns the file position of the query spectrum begins in file
 * represented as an long integer, the return value of ftell()
 *\returns -1 if failed to find the query spectrum
 */
long SpectrumCollection::binarySearchSpectrum(
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
int SpectrumCollection::matchFirstScanLine(
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

/**
 * Sets the filename of the ms2 file the spectra were parsed.
 * This function should be used only the first time the filename is set.
 * To change existing filename use set_spectrum_collection_filename.
 * Copies the value from arguement char* filename into a heap allocated memory.
 */
void SpectrumCollection::setNewFilename(
  char* filename ///< filename -in
  )
{
  int filename_length = strlen(filename) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);

  filename_ =
    strncpy(copy_filename,filename,filename_length);  
}

/**
 * sets the filename of the ms2 file the spectrum_collection was parsed
 * copies the value from arguement char* filename into a heap allocated memory
 * frees memory for the filename that is replaced
 */
void SpectrumCollection::setFilename(
  char* filename ///< filename -in
  )
{
  free(filename_);
  setNewFilename(filename);
}

/**
 * \returns the filename of the ms2 file the spectra was parsed
 * returns a char* to a heap allocated copy of the filename
 * user must free the memory
 */
char* SpectrumCollection::getFilename()
{  
  int filename_length = strlen(filename_) +1; // +\0
  char * copy_filename = 
    (char *)mymalloc(sizeof(char)*filename_length);
  return strncpy(copy_filename, filename_, filename_length);  
}

/**
 * \returns the current number of spectrum in the spectrum_collection
 */
int SpectrumCollection::getNumSpectra()
{
  return spectra_.size();
}

/**
 * \returns The current number of spectra assuming differnt charge(i.e. one spectrum with two charge states are counted as two spectra) in the spectrum_collection
 */
int SpectrumCollection::getNumChargedSpectra()
{
  return num_charged_spectra_;
}

/**
 * \returns the comments from the spectrum_collection
 * the return char* points to a newly heap allocated copy of the comments
 * user must free the new string object
 */
char* SpectrumCollection::getComment()
{
  char* comments = (char *)mycalloc(1, sizeof(char)*MAX_COMMENT);
  return strncpy(comments, comment_, MAX_COMMENT); 
}

/**
 * sets the comment of the spectrum_collection
 * copies the new_comment into a newly heap allocated copy of the comment
 */

void SpectrumCollection::setComment(
  char* new_comment ///< the new comments to be copied
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


/**
 * \returns TRUE if the spectrum_collection file has been parsed
 */
bool SpectrumCollection::getIsParsed()
{
  return is_parsed_;
}

/**
 * Takes the spectrum file name and creates a file with unique filenames.
 * The method will create one file for PSM result serializations for the 
 * target sequence and number_decoy_set number of files for decoy PSM 
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
FILE** SpectrumCollection::getPsmResultFilenames(
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
    parse_filename_path(filename_);
  
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

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
