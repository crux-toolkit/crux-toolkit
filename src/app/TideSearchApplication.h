#ifndef TIDESEARCHAPPLICATION_H
#define TIDESEARCHAPPLICATION_H

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

class TideSearchApplication : public CruxApplication {

protected:

  static bool HAS_DECOYS;

  /**
   * Free all existing mods
   */
  void cleanMods();

  void search(
    const string& spectrum_filename,
    const vector<SpectrumCollection::SpecCharge>* spec_charges,
    ActivePeptideQueue* active_peptide_queue,
    const ProteinVec& proteins,
    const vector<const pb::AuxLocation*>& locations,
    double precursor_window,
    WINDOW_TYPE_T window_type,
    double spectrum_min_mz,
    double spectrum_max_mz,
    int min_scan,
    int max_scan,
    int min_peaks,
    int search_charge,
    int top_matches,
    double highest_mz,
    OutputFiles* output_files,
    ofstream* target_file,
    ofstream* decoy_file,
    bool compute_sp,
    int aaSize, 
    double* aaFreqN,
    double* aaFreqI,
    double* aaFreqC,
    int* aaMass
  );

  void collectScoresCompiled(
    ActivePeptideQueue* active_peptide_queue,
    const Spectrum* spectrum,
    const ObservedPeakSet& observed,
    TideMatchSet::Arr2* match_arr,
    int queue_size,
    int charge
  );

  void computeWindow(
    const SpectrumCollection::SpecCharge& sc,
    WINDOW_TYPE_T window_type,
    double precursor_window,
    int max_charge,
    double* out_min,
    double* out_max,
    double* min_range,
    double* max_range
  );

  struct ScSortByMz {
    ScSortByMz(double precursor_window) { precursor_window_ = precursor_window; }
    bool operator() (const SpectrumCollection::SpecCharge x,
                     const SpectrumCollection::SpecCharge y) {
      return (x.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * x.charge <
             (y.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * y.charge;
    }
    double precursor_window_;
  };
  double bin_width_;
  double bin_offset_;

  struct InputFile {
    std::string OriginalName;
    std::string SpectrumRecords;
    bool Keep;
    InputFile(const std::string& name,
              const std::string& spectrumrecords,
              bool keep):
      OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
  };

public:

  /* 
   * This constant is the product of the original "magic number" (10000,
   * on line 4622 of search28.c) that was used to rescale the XCorr
   * score, and the integerization constant used by Benjamin Diament in
   * Tide.  In the Tide publication, that constant was reported as 10^7.
   * However, here it appears to be only 10^4.
   *
   * Note that this same constant is declared in report.cc.
   * 
   * --WSN, 10 March 2015
   */
  static const double XCORR_SCALING = 100000000.0;

  /*
   * This constant is used to put the refactored XCorr back into the
   * same range as the original XCorr score.  It is the XCorr "magic
   * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
   * tide/spectrum_preprocess2.cc).
   */
  static const double RESCALE_FACTOR = 20.0;

  bool exact_pval_search_;

  /**
   * Constructor
   */
  TideSearchApplication();

  /**
   * Destructor
   */
  ~TideSearchApplication();

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  static bool hasDecoys();

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
  virtual map<string, string> getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  virtual COMMAND_T getCommand() const;

  int calcScoreCount(
    int numelEvidenceObs,
    int* evidenceObs,
    int pepMassInt,
    int maxEvidence,
    int minEvidence,
    int maxScore,
    int minScore,
    int nAA,
    double* aaFreqN,
    double* aaFreqI,
    double* aaFreqC,
    int* aaMass,
    double* pValueScoreObs
  );

  virtual void processParams();
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
