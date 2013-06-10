#ifndef TIDESEARCHAPPLICATION_H
#define TIDESEARCHAPPLICATION_H

#include "CruxApplication.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "tide/records.h"
#include "tide/peptide.h"
#include "tide/theoretical_peak_set.h"
#include "tide/spectrum_collection.h"
#include "tide/active_peptide_queue.h"  // no include guard
#include "tide/max_mz.h"
#include "tide/report.h"

using namespace std;

typedef vector<const pb::Protein*> ProteinVec;

///////////////////////////////////////////////////////////////////////////////
// Set of peptide-spectrum matches. Supply a Reporter to Report() to see
// best few. 
class MatchSet {
public:
  typedef pair<int, int> Pair;
  typedef FixedCapacityArray<Pair> Arr;

  // Matches will be an array of pairs, (score, counter), where counter refers
  // to the index within the ActivePeptideQueue, counting from the back.  This
  // slight complication is due to the way the generated machine code fills the
  // counter in the matches buffer by decrementing the counter.
  MatchSet(Arr* matches) : matches_(matches) {}

  void Report(Reporter* reporter, int top_n, const Spectrum* spectrum,
              int charge, int spectrum_index, 
              const ActivePeptideQueue* peptides, const ProteinVec& proteins) {
    if (top_n > matches_->size())
      top_n = matches_->size();
    GetTop(top_n);
    pb::Stats stats;
    stats.set_count(matches_->size());
    /*if (FLAGS_stats) {
      double sum = 0;
      double sum_squares = 0;
      for (Arr::iterator i = matches_->begin(); i != matches_->end(); ++i) {
	double score = i->first / 100000000.0;
	sum += score;
	sum_squares += score * score;
      }
      stats.set_sum(sum);
      stats.set_sum_squares(sum_squares);
    }*/ // TODO ??
    reporter->ReportSpectrum(spectrum, charge, spectrum_index, &stats);
    for (Arr::iterator i = matches_->end() - 1; top_n; --i, --top_n) { 
      // i->second is counter as it was during scoring and corresponds to the
      // index of the peptide in the ActivePeptideQueue, counting from the
      // back. GetPeptide() retrieves the corresponding Peptide.
      const Peptide* peptide = peptides->GetPeptide(i->second);
      reporter->ReportMatch(i->first, *peptide);
    }
    reporter->WriteReport();
  }

private:
  struct less_score : public binary_function<Pair, Pair, bool> {
    // Compare scores, ignore counters.
    bool operator()(Pair x, Pair y) { return x.first < y.first; }
  };

  void GetTop(int top_n) {
    assert(top_n <= matches_->size());
    // move top n elements to end of array, with largest element last
    make_heap(matches_->begin(), matches_->end(), less_score());
    for (int i = 0; i < top_n; ++i)
      pop_heap(matches_->begin(), matches_->end() - i, less_score());
  }

  Arr* matches_;
};

class TideSearchApplication : public CruxApplication {

protected:

  void createResultsHeader(
    pb::Header* proteins_header,
    pb::Header* peptides_header,
    pb::Header* results_header,
    const string& proteins_file,
    const string& peptides_file,
    const string& spectra_file,
    double mass_window,
    int top_matches,
    const string& cmd_line
  );

  void search(
    const vector<SpectrumCollection::SpecCharge>* spec_charges,
    const ProteinVec& proteins,
    ActivePeptideQueue* active_peptide_queue,
    Reporter* reporter,
    double mass_window,
    int top_matches
  );

public:

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
  
};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
