/**
 * \file OutputFiles.h
 */
/*
 * FILE: OutputFiles.h
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Aug 24, 2009
 * PROJECT: crux
 * DESCRIPTION: A class description for handling all the various
 * output files, excluding parameter and log files.  The filenames,
 * locations and overwrite status would be taken from parameter.c.
 */
#ifndef OUTPUT_FILES_H
#define OUTPUT_FILES_H

#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include "carp.h"
#include "parameter.h"
#include "objects.h"
#include "match_collection.h"
#include "MatchFileWriter.h"

class OutputFiles{

 public:
  OutputFiles(COMMAND_T program_name);///< command printing files

  ~OutputFiles();
  void writeHeaders(int num_proteins = 0);
  void writeHeaders(const std::vector<bool>& add_this_col);
  void writeFeatureHeader(char** feature_names = NULL,
                          int num_names = 0);
  void writeFooters();
  void writeMatches(MATCH_COLLECTION_T* matches,
                    MATCH_COLLECTION_T** decoy_matches_array,
                    int num_decoys,
                    SCORER_TYPE_T rank_type = XCORR,
                    Spectrum* spectrum = NULL);
  void writeMatches(MATCH_COLLECTION_T* matches);
  void writeMatchFeatures(MATCH_T* match, 
                          double* features,
                          int num_features);

 private:
  bool createFiles(FILE*** file_array_ptr,
                   const char* output_dir,
                   const char* fileroot,
                   COMMAND_T command,
                   const char* extension,
                   bool overwrite);
  bool createFiles(MatchFileWriter*** file_array_ptr,
                   const char* output_dir,
                   const char* fileroot,
                   COMMAND_T command,
                   const char* extension);
  bool createFile(FILE** file_ptr,
                  const char* output_dir,
                  const char* filename,
                  bool overwrite);
  string makeFileName(const char* fileroot,
                      COMMAND_T command,
                      const char* target_decoy,
                      const char* extension,
                      const char* directory = NULL );
  void makeTargetDecoyList();

  void printMatchesXml(
                       MATCH_COLLECTION_T* target_matches,
                       MATCH_COLLECTION_T** decoy_matches_array,
                       Spectrum* spectrum,
                       SCORER_TYPE_T rank_type);
 

  void printMatchesTab(
    MATCH_COLLECTION_T*  target_matches, ///< from real peptides
    MATCH_COLLECTION_T** decoy_matches_array,  
                           ///< array of collections from shuffled peptides
    SCORER_TYPE_T rank_type,
    Spectrum* spectrum = NULL);

  void printMatchesSqt(
    MATCH_COLLECTION_T*  target_matches, ///< from real peptides
    MATCH_COLLECTION_T** decoy_matches_array,  
                           ///< array of collections from shuffled peptides
  Spectrum* spectrum = NULL);

  int num_files_;         ///< num files in each array
  std::string* target_decoy_list_; ///< target or decoy[-n] string of each file
  MatchFileWriter** delim_file_array_; ///< array of .txt files
  FILE** sqt_file_array_; ///< array of .sqt files
  FILE** xml_file_array_; ///< array of .xml files
  FILE*  feature_file_;   ///< file for percolator/q-ranker to write features to
  int matches_per_spec_;  ///< print this many matches per spec
  COMMAND_T command_;     ///< which crux command is writing these files
};










#endif //OUTPUT_FILES_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */






























