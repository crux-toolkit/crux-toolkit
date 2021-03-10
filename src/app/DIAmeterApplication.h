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
		if (mass_x < mass_y) { return true; }
		else if ((mass_x == mass_y) && (x.spectrum->SpectrumNumber() < y.spectrum->SpectrumNumber())) { return true; }
		else { return false; }
	}
};


class DIAmeterApplication : public CruxApplication {

 protected:
  // string output_file_name_;

  vector<InputFile> getInputFiles(const vector<string>& filepaths, int ms_level) const;
  SpectrumCollection* loadSpectra(const std::string& file);
  void loadMS1Spectra(const std::string& file, map<int, pair<double*, double*>>* ms1scan_intensity_rank_map);

  void computeWindowDIA(
		  const SpectrumCollection::SpecCharge& sc,
		  int max_charge,
		  vector<int>* negative_isotope_errors,
		  vector<double>* out_min,
		  vector<double>* out_max,
		  double* min_range,
		  double* max_range );

  double getTailorQuantile(TideMatchSet::Arr2* match_arr2);


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
 * src/./crux tide-index --peptide-list T --decoy-format peptide-reverse --missed-cleavages 2 --enzyme trypsin --max-mass 6000 --output-dir /media/ylu465/Data/proj/data/dia_search/ /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all.fasta cerevisiae_orf_trans_all
 * gdb -ex=r --args src/./crux diameter --precursor-window 10 --top-match 5 --overwrite T --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output /media/ylu465/Data/proj/data/dia_search/e01306.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all --verbosity 60 > log.txt 2> error.txt 1> output.txt
  * src/./crux tide-search --precursor-window-type mz --precursor-window 10 --top-match 5 --overwrite T --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output1 /media/ylu465/Data/proj/data/dia_search/e01306.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all
 *
 */
