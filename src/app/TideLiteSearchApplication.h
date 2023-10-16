#ifndef TIDELITESEARCHAPPLICATION_H
#define TIDELITESEARCHAPPLICATION_H

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

using namespace std;

class TideLiteSearchApplication : public CruxApplication {
 private:
 
 protected:
  static const double XCORR_SCALING;
  static const double RESCALE_FACTOR;
  static const double TAILOR_QUANTILE_TH;
  static const double TAILOR_OFFSET;

  map<pair<string, unsigned int>, bool>* spectrum_flag_;
  string output_file_name_;
  std::string remove_index_;  
  bool exact_pval_search_;
  double bin_width_;
  double bin_offset_;
  bool use_neutral_loss_peaks_;
  bool use_flanking_peaks_;

  int decoy_num_;  // Number of decoys per peptide;


  vector<InputFile> getInputFiles(const vector<string>& filepaths) const;
  static SpectrumCollection* loadSpectra(const std::string& file);
  void getPeptideIndexData(string, ProteinVec& proteins, vector<const pb::AuxLocation*>& locations, pb::Header& peptides_header);
  void createOutputFiles(std::ofstream **targe_file, std::ofstream **decoy_file);


  /**
   * Function that contains the search algorithm and performs the search
   */
  void search();

 public:

  /**
   * Constructor
   */
  TideLiteSearchApplication();

  /**
   * Destructor
   */
  ~TideLiteSearchApplication();

  /**
   * Main methods
   */
  virtual int main(int argc, char** argv);

  int main(const vector<string>& input_files);

  int main(const vector<string>& input_files, const string input_index);

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
