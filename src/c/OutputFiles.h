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
#include "MatchCollection.h"
#include "MatchFileWriter.h"
#include "PepXMLWriter.h"
#include "PinWriter.h"
#include "MzIdentMLWriter.h"

class OutputFiles{

 public:
  OutputFiles(CruxApplication* application);///< command printing files

  ~OutputFiles();
  void writeHeaders(int num_proteins = 0, bool isMixedTargetDecoy = false);
  void writeHeaders(const std::vector<bool>& add_this_col);
  void writeFeatureHeader(char** feature_names = NULL,
                          int num_names = 0);
  void writeFooters();
  void writeMatches(MatchCollection* matches,
                    std::vector<MatchCollection*>& decoy_matches_array,
                    SCORER_TYPE_T rank_type = XCORR,
                    Crux::Spectrum* spectrum = NULL);
  void writeMatches(MatchCollection* matches);
  void writeMatchFeatures(Crux::Match* match, 
                          double* features,
                          int num_features);
  void writeRankedProteins(ProteinToScore& proteinToScore,
                           MetaToRank& metaToRank,
                           ProteinToMetaProtein& proteinToMeta);
  void writeRankedPeptides(PeptideToScore& peptideToScore);
  static bool isConcat();
  static void setConcat(bool enable = true);
  static bool isProteinLevelDecoys();
  static void setProteinLevelDecoys(bool enable = true);

  bool exact_pval_search_;

 private:
  bool createFiles(FILE*** file_array_ptr,
                   const char* output_dir,
                   const char* fileroot,
                   CruxApplication* application,
                   const char* extension,
                   bool overwrite);
  bool createFiles(PepXMLWriter*** file_array_ptr,
                   const char* output_dir,
                   const char* fileroot,
                   CruxApplication* application,
                   const char* extension,
                   bool overwrite);
  bool createFiles(MatchFileWriter*** file_array_ptr,
                   const char* output_dir,
                   const char* fileroot,
                   CruxApplication* application,
                   const char* extension);

  bool createFile(MzIdentMLWriter** file_ptr,
                  const char* output_dir,
                  const char* fileroot,
                  CruxApplication* application,
                  const char* extension);

  bool createFile(FILE** file_ptr,
                  const char* output_dir,
                  const char* filename,
                  bool overwrite);

  bool createFile(
    PinWriter** pin_file_ptr,
    const char* output_dir, 
    const char* filename, 
    bool overwrite
  );
  string makeFileName(const char* fileroot,
                      CruxApplication* application,
                      const char* target_decoy,
                      const char* extension,
                      const char* directory = NULL );
  void makeTargetDecoyList();

  void printMatchesXml(
                       MatchCollection* target_matches,
                       vector<MatchCollection*>& decoy_matches_array,
                       Crux::Spectrum* spectrum,
                       SCORER_TYPE_T rank_type);
 

  void printMatchesTab(
    MatchCollection*  target_matches, ///< from real peptides
    std::vector<MatchCollection*>& decoy_matches_array,  
                           ///< array of collections from shuffled peptides
    SCORER_TYPE_T rank_type,
    Crux::Spectrum* spectrum = NULL);

  void printMatchesSqt(
    MatchCollection*  target_matches, ///< from real peptides
    std::vector<MatchCollection*>& decoy_matches_array,  
                           ///< array of collections from shuffled peptides
    Crux::Spectrum* spectrum = NULL);

  void printMatchesPin(
    MatchCollection* target_matches, ///< form real peptides 
    std::vector<MatchCollection*>& decoy_maches_array, 
                          ///< array of collection from shuffled peptides  
    Crux::Spectrum* spectrum=NULL
  );

  void printMatchesMzid(
    MatchCollection* target_matches,
    std::vector<MatchCollection*>& decoy_matches_array,
    SCORER_TYPE_T rank_type
  );

  void printMatchesMzid(
    MatchCollection* collection,
    SCORER_TYPE_T rank_type
  );

  int num_files_;         ///< num files in each array
  static bool concat_;
  static bool proteinLevelDecoys_;
  std::string* target_decoy_list_; ///< target or decoy[-n] string of each file
  MatchFileWriter** delim_file_array_; ///< array of .txt files
  FILE** sqt_file_array_; ///< array of .sqt files
  MzIdentMLWriter* mzid_file_;
 
  PepXMLWriter** xml_file_array_; ///< array of .pep.xml files
  FILE*  feature_file_;   ///< file for percolator/q-ranker to write features to 
  int matches_per_spec_;  ///< print this many matches per spec
  CruxApplication* application_;///< crux application writing these files
  PinWriter* pin_file_;///< file for inpupt of percolator
 
};










#endif //OUTPUT_FILES_H

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */






























