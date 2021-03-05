#ifndef DIAMETERAPPLICATION_H
#define DIAMETERAPPLICATION_H

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

using namespace std;


class DIAmeterApplication : public CruxApplication {

 protected:
  string output_file_name_;
  double bin_width_;
  double bin_offset_;

  static bool HAS_DECOYS;
  static bool PROTEIN_LEVEL_DECOYS;

  vector<InputFile> getInputFiles(const vector<string>& filepaths, int ms_level) const;
  SpectrumCollection* loadMS2Spectra(const std::string& file);



 public:
  static const double XCORR_SCALING;
  static const double RESCALE_FACTOR;

  /**
   * Constructor
   */
  DIAmeterApplication();

  /**
   * Destructor
   */
  ~DIAmeterApplication();


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

  virtual void processParams();


};

#endif

/*
 * src/./crux diameter --spectrum-parser pwiz --overwrite T --output-dir /media/ylu465/Data/proj/data/dia_search/diameter-output /media/ylu465/Data/proj/data/dia_search/e01306.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all-origin --verbosity 60
*  src/./crux diameter --spectrum-parser pwiz --overwrite T --output-dir /media/ylu465/Data/proj/data/dia_search/diameter-output /media/ylu465/Data/proj/data/dia_search/HYE124_TTOF5600_32fix_lgillet_L150206_001_cwt2.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all-origin --verbosity 60
 */
