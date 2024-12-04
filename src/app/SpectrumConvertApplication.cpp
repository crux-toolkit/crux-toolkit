#include <cstddef>
#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "TideIndexApplication.h"
#include "SpectrumConvertApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideMatchSet.h"
#include "tide/spectrum_collection.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> 
#include <map>
#include "tide/ActivePeptideQueue.h"
#include "residue_stats.pb.h"
#include "crux_version.h"

//here :O
//#include <iostream>
//#include <set>
/*
#include <vector>
#include <queue>
#include <sstream>
#include <cstdlib>
*/
/* This constant is the product of the original "magic number" (10000,
 * on line 4622 of search28.c) that was used to rescale the XCorr
 * score, and the integerization constant used by Benjamin Diament in
 * Tide.  In the Tide publication, that constant was reported as 10^7.
 * --WSN, 10 March 2015 */
const double SpectrumConvertApplication::XCORR_SCALING = 100000000.0;

/* This constant is used to put the refactored XCorr back into the
 * same range as the original XCorr score.  It is the XCorr "magic
 * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
 * tide/spectrum_preprocess2.cc). */
const double SpectrumConvertApplication::RESCALE_FACTOR = 20.0;

// Things done:
// -- handle terminal mod structure in the peptides internally. GH issue: #639
// -- missing charge states, 0 charege states,, override charge states handled. --> need to update doc. GH issue: #557, #607

#define CHECK(x) GOOGLE_CHECK(x)

SpectrumConvertApplication::SpectrumConvertApplication() {
/*  remove_index_ = "";
  spectrum_flag_ = NULL;
  decoy_num_ = 0;
  num_range_skipped_ = 0;
  num_precursors_skipped_ = 0;
  num_isotopes_skipped_ = 0;
  num_retained_ = 0;
  num_spectra_ = 0;
  num_spectra_searched_ = 0;
  total_candidate_peptides_ = 0;
  precursor_window_ = 0;
  spectrum_min_mz_ = 0; 
  spectrum_max_mz_ = 0;
  min_scan_ = 0;
  max_scan_ = 0;
  min_peaks_ = 0;
  min_precursor_charge_ = 0;
  max_precursor_charge_ = 0;
  out_tsv_target_  = NULL; // original tide-search output format in tab-delimited text files (txt)
  out_tsv_decoy_ = NULL;  // original tide-search output format in tab-delimited text files (txt) for the decoy psms only
  out_mztab_target_ = NULL;      // mzTAB output format
  out_mztab_decoy_ = NULL;      // mzTAB output format for the decoy psms only
  out_pin_target_ = NULL;        // pin output format for percolator
  out_pin_decoy_ = NULL;        // pin output format for percolator for the decoy psms only
  total_spectra_num_ = 0;       // The total number of spectra searched. This is counted during the spectrum conversion
*/
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {  // LOCK_TYPES are defined in model/objects.h
    locks_array_.push_back(new boost::mutex());
  }
}

SpectrumConvertApplication::~SpectrumConvertApplication() {
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    delete locks_array_[i];
  }
}

int SpectrumConvertApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}

int SpectrumConvertApplication::main(const vector<string>& input_files) {

  carp(CARP_INFO, "Running spectrum-convert...");

  num_threads_ = Params::GetInt("num-threads");
  if (num_threads_ < 1) {
    num_threads_ = boost::thread::hardware_concurrency(); // MINIMUM # = 1.
    // (Meaning just main thread) Do not make this value below 1.
  } else if (num_threads_ > 64) {
    // make sure that number of threads are reasonable, e.g. user did not specify millions of threads...
    carp(CARP_FATAL, "Requested more than 64 threads.");
  }
  carp(CARP_INFO, "Number of Threads: %d", num_threads_);


  // Convert the original file names into spectrum records if needed 
  // Update the file names in the variable inputFiles_ locally.
  // Run spectrum file convertion in parallel.
  for (vector<string>::const_iterator original_file_name = input_files.begin(); original_file_name != input_files.end(); ++original_file_name) {
    inputFiles_.push_back(SpectrumConvertApplication::InputFile(*original_file_name, *original_file_name, false));
  }
  // Launch threads to convert files
  boost::thread_group threadgroup_input_files;
  for (int t = 1; t < num_threads_; ++t) {
    boost::thread * currthread = new boost::thread(boost::bind(&SpectrumConvertApplication::getInputFiles, this, t));
    threadgroup_input_files.add_thread(currthread);
  }
  getInputFiles(0);
  // Join threads
  threadgroup_input_files.join_all();

  if (total_spectra_num_ > 0) {
    carp(CARP_INFO, "There were a total of %d spectrum conversions from %d input spectrum files.",
         total_spectra_num_, inputFiles_.size());
  }
  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);

  return 0;
}

string SpectrumConvertApplication::getOutputFileName() {
  return output_file_name_;
}

// In order to add more options, you need to add them to ./src/util/Params.cpp
vector<string> SpectrumConvertApplication::getOptions() const {
  string arr[] = {
    "fileroot",
    "num-threads",
    "output-dir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}
/*
Parameters to be removed from ./Params.cpp"
brief-output
peptide-centric-search
use-tailor-calibration
exact-p-value
charge-state
evidence-granularity
*/

string SpectrumConvertApplication::getName() const {
  return "spectrum-converter";
}

string SpectrumConvertApplication::getDescription() const {
  return
    "[[nohtml:TO BE UPDATED!!!!!! Search a collection of spectra against a sequence database, "
    "returning a collection of peptide-spectrum matches (PSMs). This is a "
    "fast search engine but requires that you first build an index with "
    "tide-index.]]"
    "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
    "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
    "</sup> algorithm, which assigns peptides to spectra by comparing the "
    "observed spectra to a catalog of theoretical spectra derived from a "
    "database of known proteins. Tide's primary advantage is its speed. Our "
    "published paper provides more detail on how Tide works. If you use Tide "
    "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
    "William Stafford Noble. <a href=\"http://dx.doi.org/10.1021/pr101196n\">"
    "&quot;Faster SEQUEST Searching for Peptide Identification from Tandem "
    "Mass Spectra&quot;</a>. <em>Journal of Proteome Research</em>. "
    "10(9):3871-9, 2011.</blockquote> "
    "<p>When <code>tide-search</code> runs, it performs "
    "several intermediate steps, as follows:</p><ol>"
    "<li>If a FASTA file was provided, convert it to an index using "
    "<code>tide-index</code>.</li>"
    "<li>Convert the given "
    "fragmentation spectra to a binary format.</li><li>Search the spectra "
    "against the database and store the results in binary format.</li><li>"
    "Convert the results to one or more requested output formats.</li></ol><p>"
    "By default, the intermediate binary files are stored in the output "
    "directory and deleted when Tide finishes execution. If you plan to search "
    "against given database more than once or search a given set of spectra "
    "more than once, then you can direct Tide to save the binary spectrum "
    "files using the <code>--store-index</code> and "
    "<code>--store-spectra</code> options. "
    "Subsequent runs of the program will go faster "
    "if provided with inputs in binary format.</p>]]";
}

vector<string> SpectrumConvertApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


vector< pair<string, string> > SpectrumConvertApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("spectrum-converter.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair("spectrum-converter.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}
bool SpectrumConvertApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T SpectrumConvertApplication::getCommand() const {
  return TIDE_SEARCH_COMMAND;  // TODO: VLAD you need to create a spectrum_convert_command
}

void SpectrumConvertApplication::processParams() {
}

void SpectrumConvertApplication::getInputFiles(int thread_id) {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  if (thread_id > inputFiles_.size())
    return;
  for (vector<SpectrumConvertApplication::InputFile>::iterator original_file_name = inputFiles_.begin()+thread_id; 
       original_file_name < inputFiles_.begin() + (inputFiles_.size()); 
       original_file_name = original_file_name + num_threads_) 
    {
    carp(CARP_DEBUG, "Start processing input files");
    bool keepSpectrumrecords = true;
    string original_name = (*original_file_name).OriginalName;
    string spectrumrecords = original_name;
    // Check if the input file is spectrum records of google protocol buffer
    pb::Header header;
    HeadedRecordReader reader(original_name, &header);
    if (header.file_type() != pb::Header::SPECTRA) {
      // converting to spectrumrecords file 

      carp(CARP_INFO, "Converting %s to spectrumrecords format", original_name.c_str());
      carp(CARP_DEBUG, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
      
      spectrumrecords = Params::GetString("store-spectra");
      keepSpectrumrecords = !spectrumrecords.empty();
      if (!keepSpectrumrecords) {
        spectrumrecords = make_file_path(FileUtils::BaseName( original_name) + ".spectrumrecords.tmp");
      } else if (inputFiles_.size() > 1) {
        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
                         "spectrum files");
      }
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      int spectra_num = 0;
      if (!SpectrumRecordWriter::convert(original_name, spectrumrecords, spectra_num)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", original_name.c_str());
      }
      locks_array_[LOCK_SPECTRUM_READING]->lock();
      total_spectra_num_ += spectra_num;
      locks_array_[LOCK_SPECTRUM_READING]->unlock();

    }
    (*original_file_name).SpectrumRecords  = spectrumrecords;
    (*original_file_name).Keep = keepSpectrumrecords;
    carp(CARP_DEBUG, "Finish converting");
  }
}

void SpectrumConvertApplication::createOutputFiles() {
  
}
