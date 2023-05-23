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
 

 protected:

  map<pair<string, unsigned int>, bool>* spectrum_flag_;
  string output_file_name_;
  std::string remove_index_;  


  vector<InputFile> getInputFiles(const vector<string>& filepaths) const;
  static SpectrumCollection* loadSpectra(const std::string& file);

  /**
   * Function that contains the search algorithm and performs the search
   */
  void search();

 public:

  // See TideLiteSearchApplication.cpp for descriptions of these two constants
  static const double XCORR_SCALING;
  static const double RESCALE_FACTOR;

  bool exact_pval_search_;

  /**
   * Constructor
   */
  TideLiteSearchApplication();

  /**
   * Destructor
   */
  ~TideLiteSearchApplication();

  /**
   * Main method
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

  virtual COMMAND_T getCommand() const;

  //void setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag);
  virtual void processParams();
  string getOutputFileName();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
