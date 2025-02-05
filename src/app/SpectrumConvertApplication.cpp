#include <cstddef>
#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
#include "SpectrumConvertApplication.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> 
#include <map>
#include "crux_version.h"
using namespace std;
#define CHECK(x) GOOGLE_CHECK(x)

SpectrumConvertApplication::SpectrumConvertApplication() {
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
  if (Params::IsDefault("spectrum-outdir")) {
    output_folder_ = Params::GetString("output-dir");
  } else {
    output_folder_ = Params::GetString("spectrum-outdir");
  }
  bool overwrite = Params::GetBool("overwrite");
  if (create_output_directory(output_folder_, overwrite)) {
    carp(CARP_FATAL, "Couldn't create output directory");
  }

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


// In order to add more options, you need to add them to ./src/util/Params.cpp
vector<string> SpectrumConvertApplication::getOptions() const {
  string arr[] = {
    "fileroot",
    "num-threads",
    "output-dir",
    "spectrum-outdir",
    "overwrite",
    "parameter-file",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

string SpectrumConvertApplication::getName() const {
  return "spectrum-converter";
}

string SpectrumConvertApplication::getDescription() const {
  return
    "[[nohtml:This command converts spectrum files into the binary spectrumrecords format]]"
    "[[html:<p>This command converts spectrum files into the binary spectrumrecords format "
    "used by the tide-search command.  Most people will not need to use this command at all, "
    "because tide-search will do the conversions automatically as needed.  However, if you plan "
    "to run multiple searches using the same input file, then you can save some time by using this "
    "command to pre-convert the spectra into spectrumrecords format.</p>]]";
}

vector<string> SpectrumConvertApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


vector< pair<string, string> > SpectrumConvertApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("<fileroot>.spectrumrecords",
   "The spectra from the input file, written in spectrumrecords format. "
   "The <fileroot> is taken from the input file."));
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
  return SPECTRUM_CONVERT_COMMAND;
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
        spectrumrecords = make_file_path(FileUtils::BaseName( original_name) + ".spectrumrecords", output_folder_);
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
