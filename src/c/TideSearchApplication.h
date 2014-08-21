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
    double* out_min,
    double* out_max
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

public:

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
  virtual string getName();

  /**
   * Returns the command description
   */
  virtual string getDescription();

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory();

  virtual COMMAND_T getCommand();

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
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
