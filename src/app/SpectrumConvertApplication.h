#ifndef TIDECONVERTAPPLICATION_H
#define TIDECONVERTAPPLICATION_H

#include "CruxApplication.h" 
#include "TideMatchSet.h"
 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "tide/theoretical_peak_set.h"
#include "tide/max_mz.h"
#include "util/MathUtil.h"
#include "tide/ActivePeptideQueue.h"
#include "TideIndexApplication.h"
#include "TideMatchSet.h"

using namespace std;



class SpectrumConvertApplication : public CruxApplication {
 private:
  struct InputFile {
    std::string OriginalName;
    std::string SpectrumRecords;
    bool Keep;
    InputFile(const std::string& name,
              const std::string& spectrumrecords,
              bool keep):
      OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
  };
 protected:

  int num_threads_;
  std::string output_folder_ = "crux-output";
  
  map<pair<string, unsigned int>, bool>* spectrum_flag_;

  int total_spectra_num_;

  vector<boost::mutex *> locks_array_;  

  //vector<SpectrumConvertApplication::InputFile> getInputFiles(const vector<string>& filepaths) const;
  void getInputFiles(int thread_id);
 
  vector<InputFile> inputFiles_;

  // sprectrum search executed in parallel threads
  
  // comparition of Spectrum data, based on neutral mass
  struct compare_spectrum{
    bool operator()(pair<pb::Spectrum, int> &spec_1, pair<pb::Spectrum, int> &spec_2){
      
      return spec_1.first.neutral_mass() > spec_2.first.neutral_mass();
    }
  };

 public:

  /**
   * Constructor
   */
  SpectrumConvertApplication();

  /**
   * Destructor
   */
  ~SpectrumConvertApplication();

  /**
   * Main methods
   */
  virtual int main(int argc, char** argv);

  int main(const vector<string>& input_files);

  /**
   * Returns the command name
   */
  virtual string getName() const;

  /**
   * Returns the command description
   */
  virtual string getDescription() const;

  /**
   * Returns the command arguments
   */
  virtual vector<string> getArgs() const;

  /**
   * Returns the command options
   */
  virtual vector<string> getOptions() const;

  /**
   * Returns the command outputs
   */
  virtual vector< pair<string, string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  /**
   * Returns the command ID 
   */
  virtual COMMAND_T getCommand() const;

  /**
   * Processes the output file names
   */
  string getOutputFileName();

  /**
   * Processes the parameters
   */
  virtual void processParams();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
