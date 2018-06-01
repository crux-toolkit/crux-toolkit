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

/**
 * Locks for multi-threading in Tide.
 */
enum _tide_search_lock {
  LOCK_RESULTS,       // Results file output
  LOCK_CASCADE,       // Only used by cascade-search on spectrum_flag (map)
  LOCK_CANDIDATES,    // Updating # of candidate peptides
  LOCK_REPORTING,     // Updating sc_index and reporting progress
  NUMBER_LOCK_TYPES   // always keep this last so the value
                      // changes as cmds are added
};

typedef enum _tide_search_lock TIDE_SEARCH_LOCK_T;

class TideSearchApplication : public CruxApplication {
private:
  //Added by Andy Lin in Feb 2016
  //function determines which mass bin a precusor mass is in
  void getMassBin (
    vector<int>& pepMassInt,
    vector<int>& pepMassIntUnique,
    ActivePeptideQueue* active_peptide_queue,
    vector<bool>* candidatePeptideStatus
  );

  //Added by Andy Lin in March 2016
  //function gets the max evidence of each mass bin(column)
  //up to mass bin of candidate precursor
  //Returns max value in curResidueEvidenceMatrix
  int getMaxColEvidence(
    const vector<vector<double> >& curResidueEvidenceMatrix,
    vector<int>& maxEvidence,
    int pepMassInt
  );

  //Added by Andy Lin in Nov 2016
  //Calculatse a residue evidence score given a
  //residue evidence matrix and a theoretical spectrum
  int calcResEvScore(
    const vector<vector<double> >& curResidueEvidenceMatrix,
    const vector<unsigned int>& intensArrayTheor,
    const vector<double>& aaMassDouble,
    Peptide* curPeptide
  );

  friend class SubtractIndexApplication;

 protected:

  struct InputFile {
    std::string OriginalName;
    std::string SpectrumRecords;
    bool Keep;
    InputFile(const std::string& name,
              const std::string& spectrumrecords,
              bool keep):
      OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
  };

  /**
  brief This variable is used with Cascade Search.
  This map contains a flag for each spectrum whether
  a spectrum has not been identified in a prior cycle (0) or not (1).
  The spectrum ID is a pair containing the ordinal number of the
  input file in the first component and a scanId-charge-state in the second component.
  The scanID and the charge state is combined into a single number
  as scanID*10 + charge_state. Charge state is required to be less than 10.
  */
  map<pair<string, unsigned int>, bool>* spectrum_flag_;
  string output_file_name_;

  static bool HAS_DECOYS;
  static bool PROTEIN_LEVEL_DECOYS;

  vector<int> getNegativeIsotopeErrors() const;
  vector<InputFile> getInputFiles(const vector<string>& filepaths) const;
  static SpectrumCollection* loadSpectra(const std::string& file);

  /**
   * Function that contains the search algorithm and performs the search
   */
  void search(void *threadarg);

  /**
    * Calls search(threadarg), and if threading, creates threads calling
    * search(threadarg)
    *
    * Call structure:
    * main -> [this function] -> search(void* threadarg)
    *                 |
    *   (if threading)|
    *                 |
    *                 -> Per Thread:
    *                           -> search(void* threadarg)
    */
  void search(
    const string& spectrum_filename,
    const vector<SpectrumCollection::SpecCharge>* spec_charges,
    vector<ActivePeptideQueue*> active_peptide_queue,
    ProteinVec& proteins,
    vector<const pb::AuxLocation*>& locations,
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
    ofstream* target_file,
    ofstream* decoy_file,
    bool compute_sp,
    int nAA,
    double* aaFreqN,
    double* aaFreqI,
    double* aaFreqC,
    int* aaMass,
    int nAARes,
    const vector<double>& dAAFreqN,
    const vector<double>& dAAFreqI,
    const vector<double>& dAAFreqC,
    const vector<double>& dAAMass,
    const pb::ModTable& mod_table,
    const pb::ModTable& nterm_mod_table,
    const pb::ModTable& cterm_mod_table,
    int numDecoys,
    vector<int>* negative_isotope_errors
  );

  void collectScoresCompiled(
    ActivePeptideQueue* active_peptide_queue,
    const Spectrum* spectrum,
    const ObservedPeakSet& observed,
    TideMatchSet::Arr2* match_arr,
    int queue_size,
    int charge
  );

  void convertResults() const;

  void computeWindow(
    const SpectrumCollection::SpecCharge& sc,
    WINDOW_TYPE_T window_type,
    double precursor_window,
    int max_charge,
    vector<int>* negative_isotope_errors,
    vector<double>* out_min,
    vector<double>* out_max,
    double* min_range,
    double* max_range
  );

  struct ScSortByMz {
    explicit ScSortByMz(double precursor_window) { precursor_window_ = precursor_window; }
    bool operator() (const SpectrumCollection::SpecCharge x,
                     const SpectrumCollection::SpecCharge y) {
      return (x.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * x.charge <
             (y.spectrum->PrecursorMZ() - MASS_PROTON - precursor_window_) * y.charge;
    }
    double precursor_window_;
  };
  double bin_width_;
  double bin_offset_;

  std::string remove_index_;

  // this map can be used to preload spectra
  // <spectrumrecords file> -> SpectrumCollection
  // the SpectrumCollection must be sorted
  std::map<std::string, SpectrumCollection*> spectra_;

 public:

  // See TideSearchApplication.cpp for descriptions of these two constants
  static const double XCORR_SCALING;
  static const double RESCALE_FACTOR;

  bool exact_pval_search_;

  /**
   * Constructor
   */
  TideSearchApplication();

  /**
   * Destructor
   */
  ~TideSearchApplication();

  unsigned int NUM_THREADS;

  /**
   * Main method
   */
  virtual int main(int argc, char** argv);

  int main(const vector<string>& input_files);

  int main(const vector<string>& input_files, const string input_index);

  static bool proteinLevelDecoys();

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

  /**
   * Struct holding necessary information for each thread to run.
   */
  struct thread_data {

    string spectrum_filename;
    const vector<SpectrumCollection::SpecCharge>* spec_charges;
    ActivePeptideQueue* active_peptide_queue;
    ProteinVec proteins;
    vector<const pb::AuxLocation*> locations;
    double precursor_window;
    WINDOW_TYPE_T window_type;
    double spectrum_min_mz;
    double spectrum_max_mz;
    int min_scan;
    int max_scan;
    int min_peaks;
    int search_charge;
    int top_matches;
    double highest_mz;
    ofstream* target_file;
    ofstream* decoy_file;
    bool compute_sp;
    int64_t thread_num;
    int64_t num_threads;
    int nAA;
    double* aaFreqN;
    double* aaFreqI;
    double* aaFreqC;
    int* aaMass;
    int nAARes;
    const vector<double>* dAAFreqN;
    const vector<double>* dAAFreqI;
    const vector<double>* dAAFreqC;
    const vector<double>* dAAMass;
    const pb::ModTable* mod_table;
    const pb::ModTable* nterm_mod_table;
    const pb::ModTable* cterm_mod_table;
    int decoysPerTarget;
    vector<boost::mutex*> locks_array;
    double bin_width;
    double bin_offset;
    bool exact_pval_search;
    map<pair<string, unsigned int>, bool>* spectrum_flag;
    int* sc_index;
    int* total_candidate_peptides;
    vector<int>* negative_isotope_errors;

    thread_data (const string& spectrum_filename_, const vector<SpectrumCollection::SpecCharge>* spec_charges_,
            ActivePeptideQueue* active_peptide_queue_, ProteinVec proteins_,
            vector<const pb::AuxLocation*> locations_, double precursor_window_,
            WINDOW_TYPE_T window_type_, double spectrum_min_mz_, double spectrum_max_mz_,
            int min_scan_, int max_scan_, int min_peaks_, int search_charge_, int top_matches_,
            double highest_mz_, ofstream* target_file_,
            ofstream* decoy_file_, bool compute_sp_, int64_t thread_num_, int64_t num_threads_, int nAA_,
            double* aaFreqN_, double* aaFreqI_, double* aaFreqC_, int* aaMass_, int nAARes_,
            const vector<double>* dAAFreqN_, const vector<double>* dAAFreqI_,
            const vector<double>* dAAFreqC_, const vector<double>* dAAMass_,
            const pb::ModTable* mod_table_, const pb::ModTable* nterm_mod_table_, const pb::ModTable* cterm_mod_table_, const int decoysPerTarget_,
            vector<boost::mutex*> locks_array_, double bin_width_, double bin_offset_, bool exact_pval_search_,
            map<pair<string, unsigned int>, bool>* spectrum_flag_, int* sc_index_, int* total_candidate_peptides_,
            vector<int>* negative_isotope_errors_) :
            spectrum_filename(spectrum_filename_), spec_charges(spec_charges_), active_peptide_queue(active_peptide_queue_),
            proteins(proteins_), locations(locations_), precursor_window(precursor_window_), window_type(window_type_),
            spectrum_min_mz(spectrum_min_mz_), spectrum_max_mz(spectrum_max_mz_), min_scan(min_scan_), max_scan(max_scan_),
            min_peaks(min_peaks_), search_charge(search_charge_), top_matches(top_matches_), highest_mz(highest_mz_),
            target_file(target_file_), decoy_file(decoy_file_), compute_sp(compute_sp_),
            thread_num(thread_num_), num_threads(num_threads_), nAA(nAA_), aaFreqN(aaFreqN_), aaFreqI(aaFreqI_), aaFreqC(aaFreqC_),
            aaMass(aaMass_), nAARes(nAARes_), dAAFreqN(dAAFreqN_), dAAFreqI(dAAFreqI_), dAAFreqC(dAAFreqC_), dAAMass(dAAMass_),
            mod_table(mod_table_), nterm_mod_table(nterm_mod_table_), cterm_mod_table(cterm_mod_table_), decoysPerTarget(decoysPerTarget_),
            locks_array(locks_array_), bin_width(bin_width_), bin_offset(bin_offset_), exact_pval_search(exact_pval_search_),
            spectrum_flag(spectrum_flag_), sc_index(sc_index_), total_candidate_peptides(total_candidate_peptides_), negative_isotope_errors(negative_isotope_errors_) {}
  };

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

  void calcResidueScoreCount (
    int nAa,
    int pepMassInt,
    vector<vector<double> >& residueEvidenceMatrix,
    vector<int>& aaMass,
    const vector<double>& aaFreqN,
    const vector<double>& aaFreqI,
    const vector<double>& aaFreqC,
    int NTermMass,
    int CTermMass,
    int minAaMass,
    int maxAaMass,
    int maxEvidence,
    int maxScore,
    vector<double>& scoreCount, //this is returned for later use
    int& scoreOffSet //this is returned for later use
  );

  double calcCombinedPval( //calculates combined p-value
    double m,
    double p,
    int numPval
  );

  int factorial(int n);

  void setSpectrumFlag(map<pair<string, unsigned int>, bool>* spectrum_flag);
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
