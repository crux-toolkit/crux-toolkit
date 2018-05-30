/**
 * \file OutputFiles.cpp
 * AUTHOR: Barbara Frewen
 * CREATE DATE: Aug 24, 2009
 * PROJECT: crux
 * \brief A class description for handling all the various
 * output files, excluding parameter and log files.
 *
 * The filenames, locations and overwrite status are taken from
 * parameter.c.
 */

#include "OutputFiles.h"
#include "util/FileUtils.h"
#include "util/Params.h"

using namespace std;
using namespace Crux;

/**
 * Default constructor for OutputFiles.  Opens all of the needed
 * files, naming them based on the values of the parameters output-dir
 * and fileroot and on the name given (search, percolator, etc.).
 * Requires that the output directory already exist. 
 */
OutputFiles::OutputFiles(CruxApplication* program_name)
: matches_per_spec_(Params::GetInt("top-match")),
  application_(program_name) {

  delim_file_array_ = NULL;
  xml_file_array_ = NULL;
  sqt_file_array_ = NULL;
  mzid_file_ = NULL;
  feature_file_ = NULL;
  pin_file_ = NULL;

  // parameters for all three file types
  bool overwrite = Params::GetBool("overwrite");
  const string output_directory = Params::GetString("output-dir");
  const string fileroot = Params::GetString("fileroot");

  int num_decoy_files = Params::GetInt("num-decoy-files");
  num_files_ = num_decoy_files + 1; // plus target file

  // TODO (BF oct-21-09): consider moving this logic to parameter.c
  COMMAND_T command = application_->getCommand();
  if (Params::GetBool("concat") || command != TIDE_SEARCH_COMMAND) {
    num_files_ = 1;
  }
  
  if (command == XLINK_SEARCH_COMMAND) {
    if (Params::GetBool("concat")) {
      num_decoy_files = 0;
      num_files_ = 1;
    } else {
      num_decoy_files = 1;
      num_files_ = 2;
    }
  }

  makeTargetDecoyList();

  carp(CARP_DEBUG, 
       "OutputFiles is opening %d files (%d decoys) in '%s' with root '%s'."
       " Overwrite: %d.", 
       num_files_, num_decoy_files, output_directory.c_str(), fileroot.c_str(), overwrite);

  // all operations create tab files
  if( Params::GetBool("txt-output") ) {
    createFiles(&delim_file_array_, 
                output_directory, 
                fileroot, 
                application_, 
                "txt");
  }

  // almost all operations create xml files
  if( command != SPECTRAL_COUNTS_COMMAND &&
      Params::GetBool("pepxml-output") ) {
    createFiles(&xml_file_array_,
                output_directory,
                fileroot,
                application_,
                "pep.xml",
                overwrite);
  }
  
  // sequest and search creates sqt files
  if( command == TIDE_SEARCH_COMMAND && Params::GetBool("sqt-output") ) {
    createFiles(&sqt_file_array_, 
                 output_directory, 
                 fileroot, 
                 application_, 
                 "sqt", 
                 overwrite);
  }

  //pin
  if ( command == TIDE_SEARCH_COMMAND && Params::GetBool("pin-output") ) {
    string filename = makeFileName(
    fileroot, 
    application_,
    NULL,// not trget and decoy file 
    "pin"
    );
    createFile(
      &pin_file_,
      output_directory, 
      filename.c_str(), 
      overwrite
    );
  }

  if( Params::GetBool("mzid-output") ) {
    createFile(&mzid_file_,
               output_directory,
               fileroot,
               application_,
               "mzid");
  }

  // only percolator and q-ranker create feature files
  if( (command == PERCOLATOR_COMMAND 
       || command == QRANKER_COMMAND)
      && Params::GetBool("feature-file") ) {
    string filename = makeFileName(fileroot, application_, 
                                   NULL, // not target or decoy
                                   "features.txt");
    createFile(&feature_file_, 
               output_directory, 
               filename.c_str(), 
               overwrite);
  }
  exact_pval_search_ = false;
}

OutputFiles::~OutputFiles() {
  for(int file_idx = 0; file_idx < num_files_; file_idx ++) {
    if( delim_file_array_ ) { delete delim_file_array_[file_idx]; }
    if( sqt_file_array_ ) { fclose(sqt_file_array_[file_idx]); }
    if( xml_file_array_ ) { xml_file_array_[file_idx]->closeFile(); }
    if(pin_file_) {pin_file_->closeFile();}
  }

  if (mzid_file_) { 
    mzid_file_->closeFile();
    delete mzid_file_;
  }

  if( feature_file_ ) { fclose(feature_file_); }

  delete [] delim_file_array_;
  delete [] sqt_file_array_;
  delete [] xml_file_array_;
  delete [] target_decoy_list_;
  delete  pin_file_;
}

/**
 * Creates an array of num_files_ strings with the target or decoy
 * tag that the file in that position should have.  The first string
 * will always be "target", the second will be "decoy" (iff num_files_
 * = 2) or "decoy-1", the third "decoy-2" and so on.
 */
void OutputFiles::makeTargetDecoyList() {
  target_decoy_list_ = new string[num_files_];
  target_decoy_list_[0] = "target";
  if( num_files_ == 2 ) {
    target_decoy_list_[1] = "decoy";
  } else {
    for(int file_idx = 1; file_idx < num_files_; file_idx++) {
      ostringstream name_builder;
      name_builder << "decoy-" << file_idx;
      target_decoy_list_[file_idx] = name_builder.str();
    }
  }
}

/**
 * \returns A string with all of the parts of the filename
 * concatenated together as
 * directory/fileroot.command-name.[target|decoy]extension.  Assumes
 * that extension includes a ".".  Either fileroot and/or target_decoy
 * may be NULL. Directory argument is optional.
 */
string OutputFiles::makeFileName(const string& fileroot,
                                 CruxApplication* application,
                                 const char* target_decoy,
                                 const char* extension,
                                 const string& directory ) {

  // get command name
  string basename_str = application->getFileStem();
  const char* basename = basename_str.c_str();

  ostringstream name_builder;
  if (!directory.empty()) {
    name_builder << directory;
    if (directory[directory.length() - 1] != '/') {
      name_builder << "/";
    }
  }
  if (!fileroot.empty()) {
    name_builder << fileroot << ".";
  }
  name_builder << basename << "." ;
  if( target_decoy != NULL && target_decoy[0] != '\0' ) {
    name_builder << target_decoy << ".";
  }
  name_builder << extension;
  string filename = name_builder.str();

  return filename;
}

/**
 * A private function for generating target and decoy files named
 * according to the given arguments.
 *
 * New files are returned via the file_array_ptr argument.  When
 * num_files > 1, exactly one target file is created and the remaining
 * are decoys.  Files are named 
 * "output-dir/fileroot.command_name.target|decoy[-n].extension".
 * Requires that the output-dir already exist and have write
 * permissions. 
 * \returns true if num_files new files are created, else false.
 */
bool OutputFiles::createFiles(FILE*** file_array_ptr,
                              const string& output_dir,
                              const string& fileroot,
                              CruxApplication* application,
                              const char* extension,
                              bool overwrite) {
  if( num_files_ == 0 ) {
    return false;
  }
  
  // allocate array
  *file_array_ptr = new FILE*[num_files_];

  // create each file
  for(int file_idx = 0; file_idx < num_files_; file_idx++ ) {
    string filename = makeFileName(fileroot, application,
                                   Params::GetBool("concat") ? NULL : target_decoy_list_[file_idx].c_str(),
                                   extension);
    createFile(&(*file_array_ptr)[file_idx], 
               output_dir, 
               filename.c_str(), 
               overwrite);

  }// next file
  
  return true;
}

/**
 * A private function for generating target and decoy pepxml files named
 * according to the given arguments.
 *
 * New files are returned via the file_array_ptr argument.  When
 * num_files > 1, exactly one target file is created and the remaining
 * are decoys.  Files are named 
 * "output-dir/fileroot.command_name.target|decoy[-n].extension".
 * Requires that the output-dir already exist and have write
 * permissions. 
 * \returns true if num_files new files are created, else false.
 */
bool OutputFiles::createFiles(PepXMLWriter*** xml_writer_array_ptr,
                              const string& output_dir,
                              const string& fileroot,
                              CruxApplication* application,
                              const char* extension,
                              bool overwrite) {
  if( num_files_ == 0 ) {
    return false;
  }
  
  // allocate array
  *xml_writer_array_ptr = new PepXMLWriter*[num_files_];

  // create each file
  for(int file_idx = 0; file_idx < num_files_; file_idx++ ) {
    string filename = makeFileName(fileroot, application,
                                   Params::GetBool("concat") ? NULL : target_decoy_list_[file_idx].c_str(),
                                   extension, output_dir);
    (*xml_writer_array_ptr)[file_idx] = new PepXMLWriter();
    (*xml_writer_array_ptr)[file_idx]->openFile(filename.c_str(), overwrite);

  }// next file
  
  return true;
}


/**
 * A private function for generating target and decoy MatchFileWriters named
 * according to the given arguments.
 *
 * MatchFileWriters are returned via the file_array_ptr argument.  When
 * num_files > 1, exactly one target file is created and the remaining
 * are decoys.  Files are named 
 * "output-dir/fileroot.command_name.target|decoy[-n].extension".
 * Requires that the output-dir already exist and have write
 * permissions. 
 * \returns true if num_files new MatchFileWriters are created, else false.
 */
bool OutputFiles::createFiles(MatchFileWriter*** file_array_ptr,
                              const string& output_dir,
                              const string& fileroot,
                              CruxApplication* application,
                              const char* extension ) {
  if( num_files_ == 0 ) {
    return false;
  }
  
  // allocate array
  *file_array_ptr = new MatchFileWriter*[num_files_];

  // create each file writer
  for(int file_idx = 0; file_idx < num_files_; file_idx++ ) {
    string filename = makeFileName(fileroot, application,
                                   Params::GetBool("concat") ? NULL : target_decoy_list_[file_idx].c_str(),
                                   extension, output_dir);
    (*file_array_ptr)[file_idx] = new MatchFileWriter(filename.c_str());
  }
  
  return true;
}


/**
 * \brief A private function for opening a file according to the given
 * arguments.
 *
 * New file is returned via the file_ptr argument.  File is named
 * output-dir/fileroot.comand_name[target_decoy].extension.  Requires that the
 * output-dir already exist and have write permissions.
 * \returns true if the file is created, else false.
 */
bool OutputFiles::createFile(FILE** file_ptr,
                             const string& output_dir,
                             const string& filename,
                             bool overwrite) {
  *file_ptr = create_file_in_path(filename, output_dir, overwrite);
  return *file_ptr != NULL;
}

bool OutputFiles::createFile(MzIdentMLWriter** file_ptr,
                             const string& output_dir,
                             const string& fileroot,
                             CruxApplication* application,
                             const char* extension) {

  string filename = makeFileName(fileroot, application, "", extension, output_dir);
  *file_ptr = new MzIdentMLWriter();
  (*file_ptr)->openFile(filename, true);
  return true;

}

/**
 * \brief A private function for opening a file according to the given
 * arguments.
 *
 * New file is returned via the file_ptr argument.  File is named
 * output-dir/fileroot.pin.  Requires that the
 * output-dir already exist and have write permissions.
 * \returns true if the file is created, else false.
 */
bool OutputFiles::createFile(
  PinWriter** pin_file_ptr,
  const string& output_dir,
  const string& filename,
  bool overwrite
) {
  *pin_file_ptr = new PinWriter();
  (*pin_file_ptr)->openFile(filename, output_dir, overwrite);  
  return pin_file_ptr != NULL;
}


/**
 * \brief Write header lines to the .txt, .sqt files, .pep.xml, and .pin
 * files.  Optional num_proteins argument for .sqt files.  Use this
 * for search commands, not post-search.
 */
void OutputFiles::writeHeaders(int num_proteins, bool isMixedTargetDecoy) {

  const char* tag = "target";

  // write headers one file at a time for tab and sqt
  for(int file_idx = 0; file_idx < num_files_; file_idx++) {
    if( delim_file_array_ ) {
      bool has_decoy = (bool)file_idx; // only first file (idx 0) is target
      if( isMixedTargetDecoy ) { // unless it is both
        has_decoy = true;
      }
        delim_file_array_[file_idx]->addColumnNames(application_, has_decoy);
        delim_file_array_[file_idx]->writeHeader();
    }

    if( sqt_file_array_ ) {
      string database = Params::GetString("protein-database");
      if (!FileUtils::Exists(database.c_str())) {
        database = Params::GetString("tide database");
      }
      MatchCollection::printSqtHeader(sqt_file_array_[file_idx],
                       tag, database, num_proteins, exact_pval_search_); 
    }
    
    if ( xml_file_array_) {
      xml_file_array_[file_idx]->writeHeader();
    }

    tag = "decoy";
  }
  //write header at a time for pin file
  if(pin_file_) {
    pin_file_->printHeader();
  }
}

/**
 * \brief Write header lines to the .txt and .pep.xml
 * files.  Use this for post-search commands, not search.
 */
void OutputFiles::writeHeaders(const vector<bool>& add_this_col) {


  // write headers one file at a time for tab and sqt
  for(int file_idx = 0; file_idx < num_files_; file_idx++) {
    if( delim_file_array_ ) {
        delim_file_array_[file_idx]->addColumnNames(application_, 
                                                    (bool)file_idx, 
                                                    add_this_col);
        delim_file_array_[file_idx]->writeHeader();
    }

    if ( xml_file_array_) {
      xml_file_array_[file_idx]->writeHeader();
    }
  }
}

/**
 * \brief Write header lines to the optional feature file.
 */
void OutputFiles::writeFeatureHeader(char** feature_names,
                                     int num_names) {
  // write feature file header
  if( feature_names && feature_file_ && num_names ) {
    fprintf(feature_file_, "scan\tlabel");
    for(int name_idx = 0; name_idx < num_names; name_idx++) {
      fprintf(feature_file_, "\t%s", feature_names[name_idx]);
    }
    fprintf(feature_file_, "\n");
  }
}

/**
 * \brief Write footer lines to .pep.xml files
 */
void OutputFiles::writeFooters() {
  if (xml_file_array_) {
    for (int file_idx = 0; file_idx < num_files_; file_idx++) {
      xml_file_array_[file_idx]->writeFooter();
    }
  }
}

/**
 * \brief Write the given matches to appropriate output files.  Limit
 * the number of matches per spectrum based on top-match parameter
 * using the ranks from rank_type.  
 */
void OutputFiles::writeMatches(
  MatchCollection*  target_matches, ///< from real peptides
  vector<MatchCollection*>& decoy_matches_array,  
                                ///< array of collections from shuffled peptides
  SCORER_TYPE_T rank_type,      ///< use ranks for this type
  Spectrum* spectrum            ///< given when all matches are to one spec
  ) {

  if (!target_matches) {
    return;  // warn?
  }

  // confirm that there are the expected number of decoy collections
  if ((int)decoy_matches_array.size() != num_files_ - 1) {
    carp(CARP_FATAL, "WriteMatches was given %d decoy collections but was expecting %d.",
         (int)decoy_matches_array.size(), num_files_ - 1);
  }

  // print to each file type
  printMatchesTab(target_matches, decoy_matches_array, rank_type, spectrum);
  printMatchesSqt(target_matches, decoy_matches_array, spectrum);
  printMatchesXml(target_matches, decoy_matches_array, spectrum, rank_type);
  printMatchesPin(target_matches, decoy_matches_array);
  printMatchesMzid(target_matches, decoy_matches_array, rank_type);
}

// already confirmed that num_files_ = num decoy collections + 1
void OutputFiles::printMatchesTab(
  MatchCollection*  target_matches, ///< from real peptides
  vector<MatchCollection*>& decoy_matches_array,  
  SCORER_TYPE_T rank_type,
  Spectrum* spectrum
) {

  carp(CARP_DETAILED_DEBUG, "Writing tab delimited results.");

  if (!delim_file_array_) {
    return;
  }

  // if a spectrum is given, use one print function
  if( spectrum ) {
    MatchCollection* cur_matches = target_matches;

    for(int file_idx = 0; file_idx < num_files_; file_idx++) {
      cur_matches->calculateDeltaCn();
      cur_matches->printTabDelimited(delim_file_array_[file_idx],
                                     matches_per_spec_, spectrum, rank_type);
      carp(CARP_DETAILED_DEBUG, "done writing file index %d", file_idx);
      if( decoy_matches_array.size() > (size_t)file_idx ) {
        cur_matches = decoy_matches_array[file_idx];
      }// else if it is NULL, num_files_ == 1 and loop will exit here
    }

  } else { // use the multi-spectra print function which assumes
           // targets and decoys are merged
    target_matches->printMultiSpectra(delim_file_array_[0],
                                (num_files_ > 1) ? delim_file_array_[1] : NULL);
  }

}

void OutputFiles::printMatchesPin(
  MatchCollection* target_matches,
  vector<MatchCollection*>& decoy_matches_array
  ) {
  if (pin_file_) {
    pin_file_->write(target_matches, decoy_matches_array, matches_per_spec_);
  }
}


void OutputFiles::printMatchesSqt(
  MatchCollection*  target_matches, ///< from real peptides
  vector<MatchCollection*>& decoy_matches_array,  
                                ///< array of collections from shuffled peptides
  Spectrum* spectrum
) {

  if( sqt_file_array_ == NULL ) {
    return;
  }

  MatchCollection* cur_matches = target_matches;

  for(int file_idx = 0; file_idx < num_files_; file_idx++) {

    cur_matches->printSqt(sqt_file_array_[file_idx],
                               matches_per_spec_,
                               spectrum);

    if( decoy_matches_array.size() > (size_t)file_idx ) {
      cur_matches = decoy_matches_array[file_idx];
    } // else if NULL, num_files_==1 and this is last loop
  }

}


void OutputFiles::printMatchesXml(
  MatchCollection*  target_matches, ///< from real peptides
  vector<MatchCollection*>& decoy_matches_array,  
                                ///< array of collections from shuffled peptides
  Spectrum* spectrum,
  SCORER_TYPE_T rank_type
  
) {
  static int index = 1;
  if (!xml_file_array_) {
    return;
  }

  MatchCollection* cur_matches = target_matches;

  for (int file_idx = 0; file_idx < num_files_; file_idx++) {
    cur_matches->printXml(xml_file_array_[file_idx], matches_per_spec_, spectrum, rank_type);
    if( decoy_matches_array.size() > (size_t)file_idx ) {
      cur_matches = decoy_matches_array[file_idx];
    } // else if NULL, num_files_==1 and this is last loop
  }
  index++;
}

void OutputFiles::printMatchesMzid(
  MatchCollection* target_matches,
  vector<MatchCollection*>& decoy_matches_array,
  SCORER_TYPE_T rank_type
  ) {
  if (!mzid_file_ || !target_matches) {
    return;
  }

  printMatchesMzid(target_matches, rank_type);

  for (size_t idx = 0; idx < decoy_matches_array.size();idx++) {
    printMatchesMzid(decoy_matches_array[idx], rank_type);
  }
}

void OutputFiles::printMatchesMzid(
  MatchCollection* collection,
  SCORER_TYPE_T rank_type
  ) {
  MatchIterator match_iter(collection, rank_type);

  while(match_iter.hasNext()) {
    Match* current_match = match_iter.next();
    if (current_match->getRank(rank_type) <= matches_per_spec_) {
      mzid_file_->addMatch(collection, current_match);
    } else {
      break;
    }
  }
}

void OutputFiles::writeMatches(
  MatchCollection*  matches ///< from multiple spectra
) {
  matches->printMultiSpectra(delim_file_array_[0], NULL /* no decoy file */);
  if (xml_file_array_) {
    matches->printMultiSpectraXml(xml_file_array_[0]);
  }
}

/**
 * \brief Print features from one match to file.
 */
void OutputFiles::writeMatchFeatures(
  Match* match, ///< match to provide scan num, decoy
  double* features,///< features for this match
  int num_features) {///< size of features array
  
  if( feature_file_ == NULL ) { return; }

  // write scan number
  fprintf(feature_file_, "%i\t",
          match->getSpectrum()->getFirstScan());

  // decoy or target peptide
  if (match->getNullPeptide() == false) {
    fprintf(feature_file_, "1\t");
  } else { 
    fprintf(feature_file_, "-1\t");
  }
  
  // print each feature, end in new-line
  for(int feature_idx = 0; feature_idx < num_features; feature_idx++) {
    if (feature_idx < num_features - 1) {
      fprintf(feature_file_, "%.4f\t", features[feature_idx]);
    } else {
      fprintf(feature_file_, "%.4f\n", features[feature_idx]);
    }
  }

}

/**
 * Print the given peptides and their scores in sorted order by score.
 */
void OutputFiles::writeRankedPeptides(const vector<pair<FLOAT_T, Peptide*> >& scoreToPeptide) {
  MatchFileWriter* file = delim_file_array_[0];
  MATCH_COLUMNS_T score_col = SIN_SCORE_COL;

  MEASURE_TYPE_T measure_type = string_to_measure_type(Params::GetString("measure"));
  switch (measure_type) {
    case MEASURE_RAW:
      score_col = RAW_SCORE_COL;
      break;
    case MEASURE_SIN:
      score_col = SIN_SCORE_COL;
      break;
    case MEASURE_NSAF:
      score_col = NSAF_SCORE_COL;
      break;
    case MEASURE_DNSAF:
      score_col = DNSAF_SCORE_COL;
      break;
    case MEASURE_EMPAI:
      score_col = EMPAI_SCORE_COL;
       break;
    default:
      carp(CARP_FATAL, "Invalid measure type!");
  }

  // print each pair
  for(vector<pair<FLOAT_T, Peptide*> >::const_iterator it = scoreToPeptide.begin();
      it != scoreToPeptide.end(); ++it) {
    Peptide* peptide = it->second;
    FLOAT_T score = it->first;
    const char* seq = peptide->getSequence();

    file->setColumnCurrentRow(SEQUENCE_COL, seq);
    file->setColumnCurrentRow(score_col, score);
    file->writeRow();
  }

}

void OutputFiles::pinSetEnabledStatus(const string& name, bool enabled) {
  if (pin_file_) {
    pin_file_->setEnabledStatus(name, enabled);
  }
}

/**
 * Print all of the proteins and their associated scores in sorted
 * order by score. If there is parsimony information, also print the
 * parsimony rank.
 */
void OutputFiles::writeRankedProteins(const vector<boost::tuple<FLOAT_T, Protein*, int> >& proteins,
                                      bool isParsimony) {
  MatchFileWriter* file = delim_file_array_[0];
  MATCH_COLUMNS_T score_col = SIN_SCORE_COL;

  switch (string_to_measure_type(Params::GetString("measure"))) {
    case MEASURE_RAW:
      score_col = RAW_SCORE_COL;
      break;
    case MEASURE_SIN:
      score_col = SIN_SCORE_COL;
      break;
    case MEASURE_NSAF:
      score_col = NSAF_SCORE_COL;
      break;
    case MEASURE_DNSAF:
      score_col = DNSAF_SCORE_COL;
      break;
    case MEASURE_EMPAI:
      score_col = EMPAI_SCORE_COL;
       break;
    default:
      carp(CARP_FATAL, "Invalid measure type!");
  }

  // print each protein
  for(vector<boost::tuple<FLOAT_T, Protein*, int> >::const_iterator it = proteins.begin();
      it != proteins.end(); ++it) {
    if (isParsimony && it->get<2>() < 1) {
      // Don't print out -1 ranked proteins
      carp(CARP_DEBUG, "-1 Parsimony Rank for protein detected:%s",
	   it->get<1>()->getIdPointer().c_str());
      continue;
    }
    file->setColumnCurrentRow(score_col, it->get<0>());
    file->setColumnCurrentRow(PROTEIN_ID_COL, it->get<1>()->getIdPointer());
    if (isParsimony) {
      file->setColumnCurrentRow(PARSIMONY_RANK_COL, it->get<2>());
    }
    file->writeRow();
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
