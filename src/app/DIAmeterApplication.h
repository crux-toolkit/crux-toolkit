/**
 * DIAmeterApplication.h
 * DATE: June 15, 2021
 * AUTHOR: Yang Lu
 * DESCRIPTION: The DIAmeter module for DIA search.
 **************************************************************************/

#ifndef DIAMETERAPPLICATION_H
#define DIAMETERAPPLICATION_H

#include "CruxApplication.h"
#include "TideMatchSet.h"
#include "TideSearchApplication.h"

#include <iostream> 
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "tide/theoretical_peak_set.h"
#include "tide/max_mz.h"

using namespace std;
struct InputFile {
  std::string OriginalName;
  std::string SpectrumRecords;
  bool Keep;
  InputFile(const std::string& name,
            const std::string& spectrumrecords,
            bool keep):
    OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
};
// It is modified from TideSearchApplication ScSortByMz with following differences:
// By default the mz tolerane would be half of the isolation window
// In case the isolation window is 0 which should cause a fatal error elsewhere
// We first want the spectra is sorted by neutral mass like usual.
// Then if the neutral mass is equal, sort by the MS2 scan
struct ScSortByMzDIA {
    bool operator() (const SpectrumCollection::SpecCharge x, const SpectrumCollection::SpecCharge y) {
        double precursor_window_x = fabs(x.spectrum->IsoWindowUpperMZ()-x.spectrum->IsoWindowLowerMZ()) / 2;
        double precursor_window_y = fabs(y.spectrum->IsoWindowUpperMZ()-y.spectrum->IsoWindowLowerMZ()) / 2;
        // carp(CARP_DETAILED_DEBUG, "precursor_window_x:%f \t precursor_window_y:%f", precursor_window_x, precursor_window_y);
        //return (x.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_x) * x.charge < (y.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_y) * y.charge;

        double mass_x = (x.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_x) * x.charge;
        double mass_y = (y.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_y) * y.charge;
        if (mass_x < mass_y) {
          return true;
        } else if ((mass_x == mass_y) && (x.spectrum->SpectrumNumber() < y.spectrum->SpectrumNumber())) {
          return true;
        } else {
          return false;
        }
    }
};

class DIAmeterApplication : public CruxApplication {

 protected:
  // string output_file_name_;
  double avg_noise_intensity_logrank_, avg_ms1_intercept_, avg_isowin_width_;
  int scan_gap_, max_ms1scan_;

  std::string remove_index_, output_pin_, output_percolator_;

  vector<InputFile> getInputFiles(const vector<string>& filepaths, int ms_level) const;

  SpectrumCollection* loadSpectra(const std::string& file);

  void loadMS1Spectra(const std::string& file,
          map<int, boost::tuple<double*, double*, double*, int>>* ms1scan_mz_intensity_rank_map,
          map<int, boost::tuple<double, double>>* ms1scan_slope_intercept_map
  );

  void buildSpectraIndexFromIsoWindow(vector<SpectrumCollection::SpecCharge>* spec_charge_chunk, map<int, boost::tuple<double*, double*, int>>* ms2scan_mz_intensity_map);

  void reportDIA(
    ofstream* output_file,  // output file to write to
    const string& spectrum_filename, // name of spectrum file
    const SpectrumCollection::SpecCharge& sc, // spectrum and charge for matches
    ActivePeptideQueue* peptides, // peptide queue
    const ProteinVec& proteins, // proteins corresponding with peptides
    TideMatchSet& matches, // object to manage PSMs
    ObservedPeakSet* observed,
    map<int, boost::tuple<double*, double*, double*, int>>* ms1scan_mz_intensity_rank_map,
    map<int, boost::tuple<double, double>>* ms1scan_slope_intercept_map,
    map<int, boost::tuple<double*, double*, int>>* ms2scan_mz_intensity_map,
    map<string, double>* peptide_predrt_map
  );

  void computePrecIntRank(
      TideMatchSet::PSMScores& vec,
      ActivePeptideQueue* peptides,
      const double* mz_arr,
    const double* intensity_arr,
      const double* intensity_rank_arr,
    boost::tuple<double, double> slope_intercept_tp,
    int peak_num,
      map<TideMatchSet::PSMScores::iterator, boost::tuple<double, double, double>>* intensity_map,
      map<TideMatchSet::PSMScores::iterator, boost::tuple<double, double, double>>* logrank_map,
      int charge
  );

  void computePrecFragCoelute(
    TideMatchSet::PSMScores& vec,
    ActivePeptideQueue* peptides,
    vector<boost::tuple<double*, double*, int, double*, double*, int>>* mz_intensity_arrs_vector,
      map<TideMatchSet::PSMScores::iterator, boost::tuple<double, double, double>>* coelute_map,
      int charge
  );

  void computeMS2Pval(
    TideMatchSet::PSMScores& vec,
    ActivePeptideQueue* peptides,
    ObservedPeakSet* observed,
    map<TideMatchSet::PSMScores::iterator, boost::tuple<double, double>>* ms2pval_map
  );

  void computeWindowDIA(
    const SpectrumCollection::SpecCharge& sc,
    vector<int>* negative_isotope_errors,
    vector<double>* out_min,
    vector<double>* out_max,
    double* min_range,
    double* max_range
  );

  void getPeptidePredRTMapping(map<string, double>* peptide_predrt_map, int percent_bins = 200);

  double closestPPMValue(
    const double* mz_arr,
    const double* intensity_arr,
    int peak_num,
    double query_mz,
    int ppm_tol,
    double intensity_default,
    bool large_better
  );

  string getCoeffTag();

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

