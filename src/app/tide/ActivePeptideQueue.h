#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"
#include "IonInvertedIndex.h"

#ifndef ACTIVE_PEPTIDE_QUEUE_H
#define ACTIVE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

// Constants for E-value calculation from score distribution (all scores, target + decoy)
static const int EVALUE_HISTOGRAM_SIZE = 2000;  // max xCorr * 10; bin count for E-value histogram
static const int MIN_SCORES_FOR_REGRESSION = 2;  // minimum scores required to fit a regression model for E-values
static const double EVALUE_UPPER_BOUND = 999.0;  // maximum reported E-value

class ActivePeptideQueue {
 public:
  ActivePeptideQueue(RecordReader* reader,
                     const vector<const pb::Protein*>& proteins,
                     vector<const pb::AuxLocation*>* locations=NULL,
                     bool dia_mode = false);

  ~ActivePeptideQueue();

  int SetActiveRange(double min_range, double max_range, double max_exp_peak_mz);
  inline size_t size() const { return queue_.size(); }

  Peptide* GetPeptide(int index) {
    size_t sz = end_ - begin_;
    if (index >= static_cast<int>(sz)) index = static_cast<int>(sz) - 1;
    if (index < 0) index = 0;
    return *(begin_ + index);
  }

  SCORE_FUNCTION_T curScoreFunction_;

  // Peptide counting statistics
  int total_peptides_loaded_;     // total peptides currently in the queue
  int candidate_peptide_count_;   // peptides within the precursor mass tolerance window
  int candidate_target_count_;    // candidate peptides that are targets
  int candidate_decoy_count_;     // candidate peptides that are decoys

  // --- Tailor score calibration histogram ---
  // Collects scores from ALL scored peptides (target + decoy) in both XCorr and
  // HyperScore modes. Used for:
  //   1. Tailor quantile calibration (GetTailorQuantile) — applied to both XCorr
  //      and HyperScore PSMs.
  //   2. Poisson E-value for HyperScore — GetTailorMeanMatches() provides the
  //      lambda parameter (average ion matches per peptide).
  vector<int> tailor_histogram_;

  // Histogram discretization: bin = round((score + tailor_histogram_offset_) * histogram_bin_scale_)
  int histogram_bin_scale_;          // scale factor converting real-valued scores to integer bins
  int tailor_histogram_offset_;      // offset to ensure non-negative bin indices (needed for negative XCorr scores)
  int tailor_histogram_max_score_;   // maximum score the histogram can represent (before resizing)
  int tailor_histogram_max_bin_;     // highest bin index that has received a score
  int tailor_score_count_;           // number of scores added to the Tailor histogram
  int tailor_total_matches_;         // total ion matches across all scored peptides

  // Score → bin conversion:
  //   int bin = round((score + tailor_histogram_offset_) * histogram_bin_scale_);
  //   ++tailor_histogram_[bin];
  // Bin → score conversion:
  //   score = (double)bin / histogram_bin_scale_ - tailor_histogram_offset_

  const double TAILOR_QUANTILE_TH = 0.01;
  const double TAILOR_OFFSET = 5.0;

  // Tailor calibration: returns the score at the 99th percentile of the Tailor histogram
  double GetTailorQuantile();
  // Returns the average ion matches per scored peptide; used as Poisson
  // lambda for the HyperScore Poisson-based E-value calculation.
  double GetTailorMeanMatches();

  void ResetTailorHistogram();
  void AddToTailorHistogram(double score, int match_cnt = 0);

  // --- E-value calculation from score distribution (all scores, target + decoy) ---
  // The E-value histogram collects ALL peptide scores (not just decoys).
  // For XCorr:  regression on the cumulative score distribution → ComputeEValue(xcorr)
  // For HyperScore: regression on the survival function S(x) = P(score >= x) → ComputeEValue(hyper)
  // HyperScore also supports a Poisson-based E-value using GetTailorMeanMatches() as lambda.
  void ResetEValueHistogram();   // clear E-value histogram before scoring a spectrum
  void AddToEValueHistogram(double score);  // record a score for E-value computation
  void FitEValueRegression();    // fit regression model after all scores are collected
  double ComputeEValue(double score) const;  // E-value from fitted regression
  int ScoreToBin(double score) const;  // map a score to an E-value histogram bin
  void LinearRegression(const int* histogram, double* slope, double* intercept,
                        int* max_bin, int* start_bin, int* end_bin);
  bool LinearRegressionHyperScore(const int* histogram, int hist_size,
                                  double* slope, double* intercept, double* rsq);
  // class PeptideWrapper {  // Peptide objects split into hot and cold data in order to reduce cache miss ratio. 
  // hot data:
  //   double mass;
  //   bool isDecoy;
  //   vector<unsigned int> peaks_0;   // Single charged b-y ions, in case of exact p-value, this contains only the b-ions
  //   vector<unsigned int> peaks_1;   // Double charged b-y ions
  // cold data:
  // Peptide* peptide;
  // }

  deque<Peptide*> queue_;
  deque<Peptide*>::const_iterator begin_, end_;

  int min_candidates_;       // minimum peptides to buffer for Tailor calibration
  bool dia_mode_;            // DIA (data-independent acquisition) mode flag

  IonInvertedIndex ion_inverted_index_;

 private:

  void ComputeTheoreticalPeaksBack();

  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_;
  vector<const pb::AuxLocation*>* locations_;

  TheoreticalPeakSetBYSparse theoretical_peak_set_;
  pb::Peptide current_pb_peptide_;

  // --- E-value regression state (per-spectrum) ---
  int evalue_histogram_[EVALUE_HISTOGRAM_SIZE];  // histogram of ALL scores for E-value; index i → score = i / 10.0
  int evalue_score_count_;                       // number of scores added for current spectrum

  // Linear regression parameters fitted on the score distribution
  double evalue_slope_;
  double evalue_intercept_;
  int evalue_reg_start_bin_;   // lower bound of regression range (in bins)
  int evalue_reg_end_bin_;     // upper bound of regression range (in bins)
  int evalue_max_bin_;         // highest bin with a non-zero count
};

#endif
