/**
 * \file output-files.cpp
 */
/*
 * FILE: output-files.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Aug 24, 2009
 * PROJECT: crux
 * DESCRIPTION: A class description for handling all the various
 * output files, excluding parameter and log files.  The filenames,
 * locations and overwrite status would be taken from parameter.c.
 */

#include "output-files.h"
using namespace std;

/**
 * Default constructor for OutputFiles.  Opens all of the needed
 * files, naming them based on the values of the parameters output-dir
 * and fileroot and on the name given (search, percolator, etc.).
 * Requires that the output directory already exist. 
 */
OutputFiles::OutputFiles(const char* program_name)
: matches_per_spec_(get_int_parameter("top-match"))
{

  psm_file_array_ = NULL;
  tab_file_array_ = NULL;
  sqt_file_array_ = NULL;

  // parameters for all three file types
  BOOLEAN_T overwrite = get_boolean_parameter("overwrite");
  const char* output_directory = get_string_parameter_pointer("output-dir");
  const char* fileroot = get_string_parameter_pointer("fileroot");
  if( strcmp(fileroot, "__NULL_STR") == 0 ){
    fileroot = NULL;
  }
  int num_decoy_files = get_int_parameter("num-decoy-files");
  num_files_ = num_decoy_files + 1; // plus target file

  carp(CARP_DEBUG, 
       "OutputFiles is opening %d files (%d decoys) in '%s' with root '%s'."
       " Overwrite: %d.", 
       num_files_, num_decoy_files, output_directory, fileroot, overwrite);

  // all operations create tab files
  createFiles(&tab_file_array_, 
              output_directory, 
              fileroot, 
              program_name, 
              "txt", 
              overwrite); 


  // search operations create .csm files
  if( strcmp(program_name, "search") == 0 ||
      strcmp(program_name, "sequest") == 0 ){
    createFiles(&psm_file_array_, 
                 output_directory, 
                 fileroot, 
                 program_name, 
                 "csm", 
                 overwrite);
  }

  // only sequest creates sqt files
  //if( strcmp(program_name, "sequest") == 0 ){
  // TEMPORARY: also write sqt for search-for-matches so that smoke
  //tests will pass
  if( strcmp(program_name, "sequest") == 0 || strcmp(program_name, "search") == 0){
    createFiles(&sqt_file_array_, 
                 output_directory, 
                 fileroot, 
                 program_name, 
                 "sqt", 
                 overwrite);
  }
}

OutputFiles::~OutputFiles(){
  delete [] psm_file_array_;
  delete [] tab_file_array_;
  delete [] sqt_file_array_;
}

/**
 * A private function for generating target and decoy files named
 * according to the given arguments.
 *
 * New files are returned via the file_array_ptr argument.  When
 * num_files > 1, exactly one target file is created and the remaining
 * are decoys.  Files are named 
 * "output-dir/fileroot.command_name.target|decoy[n].extension".
 * Requires that the output-dir already exist and have write
 * permissions. 
 * \returns TRUE if num_files new files are created, else FALSE.
 */
BOOLEAN_T OutputFiles::createFiles(FILE*** file_array_ptr,
                                   const char* output_dir,
                                   const char* fileroot,
                                   const char* command_name,
                                   const char* extension,
                                   BOOLEAN_T overwrite){
  if( num_files_ == 0 ){
    return FALSE;
  }
  
  // allocate array
  *file_array_ptr = (FILE**)mycalloc(num_files_, sizeof(FILE*));

  // determine the name for each of the files: target, decoy-1, etc.
  string* target_decoy = new string[num_files_];
  target_decoy[0] = "target";
  if( num_files_ == 2 ){
    target_decoy[1] = "decoy";
  }else{
    for(int file_idx = 1; file_idx < num_files_; file_idx++){
      ostringstream name_builder;
      name_builder << "decoy-" << file_idx;
      target_decoy[file_idx] = name_builder.str();
    }
  }
  
  // create each file
  for(int file_idx = 0; file_idx < num_files_; file_idx++ ){
    // concatinate the pieces of the name
    ostringstream name_builder;
    if( fileroot ){
      name_builder << fileroot << "." ;
    }
    name_builder << command_name << "."
                 << target_decoy[file_idx] << "." 
                 << extension;
    string filename = name_builder.str();

    // open the file (it checks for success)
    (*file_array_ptr)[file_idx] = create_file_in_path(filename.c_str(),
                                                      output_dir,
                                                      overwrite);
  }// next file

  delete [] target_decoy;

  return TRUE;
}

void OutputFiles::writeHeaders(int num_proteins){

  const char* tag = "target";

  // write headers one file at a time for tab and sqt
  for(int file_idx = 0; file_idx < num_files_; file_idx++){
    if( tab_file_array_ ){
      print_tab_header(tab_file_array_[file_idx]);
    }

    if( sqt_file_array_ ){
      print_sqt_header(sqt_file_array_[file_idx],
                       tag,
                       num_proteins, FALSE); // not post search
    }
    tag = "decoy";
  }
  // write headers for all three psm files at once
  if( psm_file_array_ ){
    serialize_headers(psm_file_array_);
  }

}
/**
 * \brief Write the given matches to appropriate output files.  Limit
 * the number of matches per spectrum based on top-match parameter
 * using the ranks from rank_type.  
 */
// TODO: ensure this causes no changes to the match collections by
// making them const
void OutputFiles::writeMatches(
  MATCH_COLLECTION_T*  target_matches, ///< from real peptides
  MATCH_COLLECTION_T** decoy_matches_array,  
                                ///< array of collections from shuffled peptides
  int num_decoy_collections,    ///< num collections in array
  SCORER_TYPE_T rank_type,           ///< use ranks for this type
  SPECTRUM_T* spectrum     ///< given when all matches are to one spec
  ){

  if( target_matches == NULL ){
    return;  // warn?
  }

  // confirm that there are the expected number of decoy collections
  if( num_decoy_collections != num_files_ - 1){
    carp(CARP_FATAL, 
         "WriteMatches was given %d decoy collections but was expecting %d.",
         num_decoy_collections, num_files_ - 1);
  }

  // print to each file type
  printMatchesTab(target_matches, decoy_matches_array, rank_type, spectrum);
  
  printMatchesPsm(target_matches, decoy_matches_array);
  
  printMatchesSqt(target_matches, decoy_matches_array, spectrum);

}

// already confirmed that num_files_ = num decoy collections + 1
void OutputFiles::printMatchesTab(
  MATCH_COLLECTION_T*  target_matches, ///< from real peptides
  MATCH_COLLECTION_T** decoy_matches_array,  
  SCORER_TYPE_T rank_type,
  SPECTRUM_T* spectrum
){

  carp(CARP_DETAILED_DEBUG, "Writing tab delimited results.");

  if( tab_file_array_ == NULL ){
    return;
  }

  // if a spectrum is given, use one print function
  if( spectrum ){
    MATCH_COLLECTION_T* cur_matches = target_matches;

    for(int file_idx = 0; file_idx < num_files_; file_idx++){

      print_match_collection_tab_delimited(tab_file_array_[file_idx],
                                           matches_per_spec_,
                                           cur_matches,
                                           spectrum,
                                           rank_type);

      if( decoy_matches_array ){
        cur_matches = decoy_matches_array[file_idx];
      }// else if it is NULL, num_files_ == 1 and loop will exit here
    }

  } else { // use the multi-spectra print function which assumes
           // targets and decoys are merged
    print_matches_multi_spectra(target_matches,
                                tab_file_array_[0],
                                (num_files_ > 1) ? tab_file_array_[1] : NULL);
  }

}

void OutputFiles::printMatchesPsm(
  MATCH_COLLECTION_T*  target_matches, ///< from real peptides
  MATCH_COLLECTION_T** decoy_matches_array  
                                ///< array of collections from shuffled peptides
){
  if( psm_file_array_ == NULL ){
    return;
  }

  carp(CARP_DETAILED_DEBUG, "Writing csm results.");

  MATCH_COLLECTION_T* cur_matches = target_matches;
  
  for(int file_idx = 0; file_idx < num_files_; file_idx++){
    
    serialize_psm_features(cur_matches,
                           psm_file_array_[file_idx],
                           matches_per_spec_,
                           SP, XCORR);
    
    if( decoy_matches_array ){
      cur_matches = decoy_matches_array[file_idx];
    } // else if NULL, num_files_==1 and this is last loop
  }
  
}

void OutputFiles::printMatchesSqt(
  MATCH_COLLECTION_T*  target_matches, ///< from real peptides
  MATCH_COLLECTION_T** decoy_matches_array,  
                                ///< array of collections from shuffled peptides
  SPECTRUM_T* spectrum
){

  if( sqt_file_array_ == NULL ){
    return;
  }

  MATCH_COLLECTION_T* cur_matches = target_matches;

  for(int file_idx = 0; file_idx < num_files_; file_idx++){

    print_match_collection_sqt(sqt_file_array_[file_idx],
                               matches_per_spec_,
                               cur_matches,
                               spectrum,
                               SP, XCORR);

    if( decoy_matches_array ){
      cur_matches = decoy_matches_array[file_idx];
    } // else if NULL, num_files_==1 and this is last loop
  }

}

void OutputFiles::updateHeaders(int spectrum_count){

  if( psm_file_array_ == NULL ){
    return;
  }

  for(int file_idx = 0; file_idx < num_files_; file_idx++){
    serialize_total_number_of_spectra(spectrum_count,
                                      psm_file_array_[file_idx]);
  }
}



















/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
