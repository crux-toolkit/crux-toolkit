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
  double avg_noise_intensity_logrank_;
  int scan_gap_;
  int max_ms1scan_;
  int max_ms1_mzbin_, max_ms2_mzbin_;


  vector<InputFile> getInputFiles(const vector<string>& filepaths, int ms_level) const;

  SpectrumCollection* loadSpectra(const std::string& file);

  void loadMS1Spectra(const std::string& file, map<int, pair<double*, double*>>* ms1scan_intensity_rank_map);

  void buildSpectraIndexFromIsoWindow(vector<SpectrumCollection::SpecCharge>* spec_charge_chunk, map<int, double*>* ms2scan_intensity_map);

  void reportDIA(
	ofstream* output_file,  ///< output file to write to
	const string& spectrum_filename, ///< name of spectrum file
	const SpectrumCollection::SpecCharge& sc, ///< spectrum and charge for matches
	const ActivePeptideQueue* peptides, ///< peptide queue
	const ProteinVec& proteins, ///< proteins corresponding with peptides
	const vector<const pb::AuxLocation*>& locations,  ///< auxiliary locations
	TideMatchSet* matches, ///< object to manage PSMs
	ObservedPeakSet* observed,
	map<int, pair<double*, double*>>* ms1scan_intensity_rank_map,
	map<int, double*>* ms2scan_intensity_map,
	map<string, double>* peptide_predrt_map
  );

  void computePrecIntRank(
	const vector<TideMatchSet::Arr::iterator>& vec,
	const ActivePeptideQueue* peptides,
	const double* intensity_rank_arr,
	map<TideMatchSet::Arr::iterator, boost::tuple<double, double, double>>* intensity_map,
	int charge
  );

  void computeMS2Pval(
	const vector<TideMatchSet::Arr::iterator>& vec,
	const ActivePeptideQueue* peptides,
	ObservedPeakSet* observed,
	map<TideMatchSet::Arr::iterator, double>* ms2pval_map,
	bool dynamic_filter = true
  );

  void computePrecFragCoelute(
  	const vector<TideMatchSet::Arr::iterator>& vec,
  	const ActivePeptideQueue* peptides,
	vector<pair<double*, double*>>* intensity_arrs_vector,
	map<TideMatchSet::Arr::iterator, boost::tuple<double, double, double>>* coelute_map,
	int charge
  );


  void computeWindowDIA(
		  const SpectrumCollection::SpecCharge& sc,
		  int max_charge,
		  vector<int>* negative_isotope_errors,
		  vector<double>* out_min,
		  vector<double>* out_max,
		  double* min_range,
		  double* max_range
  );

  double getTailorQuantile(TideMatchSet::Arr2* match_arr2);

  void getPeptidePredRTMapping(map<string, double>* peptide_predrt_map, int percent_bins=200);


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
 * TODO remove temporary commands finally
 * src/./crux tide-index --peptide-list T --decoy-format peptide-reverse --missed-cleavages 2 --enzyme trypsin --max-mass 6000 --output-dir /media/ylu465/Data/proj/data/dia_search/ /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all.fasta cerevisiae_orf_trans_all
 * gdb -ex=r --args src/./crux diameter --precursor-window 10 --precursor-window-type mz --top-match 5 --use-tailor-calibration T --concat T --overwrite T --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output /media/ylu465/Data/proj/data/dia_search/e01306.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all --predrt-files /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all/deeprt.peptides.target.txt,/media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all/deeprt.peptides.decoy.txt --verbosity 60 > log.txt 2> error.txt 1> output.txt
 * src/./crux tide-search --precursor-window 10 --precursor-window-type mz --use-tailor-calibration T --top-match 5 --concat T --overwrite T --num-threads 1 --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output /media/ylu465/Data/proj/data/dia_search/e01306.mzXML /media/ylu465/Data/proj/data/dia_search/cerevisiae_orf_trans_all
 *
 * src/./crux make-pin --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output --output-file diameter-search.filtered.pin --overwrite T /media/ylu465/Data/proj/data/dia_search/crux-output/diameter-search.filtered.txt
 * src/./crux percolator --tdc F --output-weights T --overwrite T --unitnorm T --pepxml-output T --output-dir /media/ylu465/Data/proj/data/dia_search/crux-output/percolator_prec_1.00_frag_1.00_rt_1.00_elu_1.00 /media/ylu465/Data/proj/data/dia_search/crux-output/diameter-search.filtered_prec_1.00_frag_1.00_rt_1.00_elu_1.00.txt.pin
 *
 */
