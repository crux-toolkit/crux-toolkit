// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
#include <algorithm>
#include <deque>
#include <gflags/gflags.h>
#include "records.h"
#include "peptides.pb.h"
#include "peptide.h"
#include "ActivePeptideQueue.h"
#include "records_to_vector-inl.h"
#include "theoretical_peak_set.h"
#include "compiler.h"
#include "app/TideMatchSet.h"
#include <map> 
#define CHECK(x) GOOGLE_CHECK((x))
// #define HISTOGRAM_SIZE  2000   // scores 0 .. 1999
#define MIN_TAIL_POINTS    3   // minimum points required for a valid regression
#define GAP_THRESHOLD      5   // consecutive zero bins that terminate the tail
#define SURV_MAX        0.10   // upper S(x) bound: skip bulk where S > 10%


ActivePeptideQueue::ActivePeptideQueue(RecordReader* reader,
                                       const vector<const pb::Protein*>& proteins,
                                       vector<const pb::AuxLocation*>* locations,
                                       bool dia_mode)
  : reader_(reader),
    proteins_(proteins),
    theoretical_peak_set_(1000),   // probably overkill, but no harm
    locations_(locations),
    dia_mode_(dia_mode) {
  CHECK(reader_->OK());
  min_candidates_ = 1000;
  total_peptides_loaded_ = 0;
  candidate_peptide_count_ = 0;
  candidate_target_count_ = 0;
  candidate_decoy_count_ = 0;


  tailor_histogram_offset_ = 10;   // Sufficient offset for XCorr (which can be negative); irrelevant for HyperScore
  histogram_bin_scale_ = 100;   // TODO: may be fine-tuned experimentally.
  tailor_histogram_max_score_ = 100;
  tailor_histogram_max_bin_ = 0;
  tailor_score_count_ = 1; // ensure the score histogram vector is allocated
  tailor_histogram_.resize((tailor_histogram_max_score_ + tailor_histogram_offset_) * histogram_bin_scale_);
  ResetTailorHistogram();
  curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function"));

  // E-value fields initialization
  memset(evalue_histogram_, 0, sizeof(evalue_histogram_));
  evalue_score_count_ = 0;
  evalue_slope_ = 0.0;
  evalue_intercept_ = 0.0;
  evalue_reg_start_bin_ = 0;
  evalue_reg_end_bin_ = 0;
  evalue_max_bin_ = 0;
}

void ActivePeptideQueue::AddToTailorHistogram(double score, int match_cnt) {
  int bin = round((score + tailor_histogram_offset_) * histogram_bin_scale_);
  tailor_histogram_[bin] += 1;
  ++tailor_score_count_;
  tailor_total_matches_ += match_cnt;
  if (tailor_histogram_max_bin_ < bin)
    tailor_histogram_max_bin_ = bin;
}

void ActivePeptideQueue::ResetTailorHistogram() {
  tailor_total_matches_ = 0;
  tailor_histogram_max_bin_ = 0;
  if (tailor_score_count_ == 0)
    return;
  //memset(tailor_histogram_.data(), 0, tailor_histogram_.capacity() * sizeof(int));
  std::fill(tailor_histogram_.begin(), tailor_histogram_.end(), 0); // reset the values of the score histogram
  tailor_score_count_ = 0;
}

double ActivePeptideQueue::GetTailorQuantile() {
  int bin = tailor_histogram_max_bin_ + 1;
  int count = 0;
  double quantile_score = 1.0;

  int quantile_pos = (int)(TAILOR_QUANTILE_TH * (double)tailor_score_count_ + 0.5) - 1; // zero indexed
  if (quantile_pos < 2)
    quantile_pos = 2;  // the third element
  if (quantile_pos >= tailor_score_count_)
    quantile_pos = tailor_score_count_ - 1; // the last element

  while (count < quantile_pos) {
    count += tailor_histogram_[bin--];
  }
  ++bin;
  return (double)bin / histogram_bin_scale_ - tailor_histogram_offset_;
}

double ActivePeptideQueue::GetTailorMeanMatches() {
  return (double)tailor_total_matches_ / (double)tailor_score_count_;
}

ActivePeptideQueue::~ActivePeptideQueue() {
}

// Compute the theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueue::ComputeTheoreticalPeaksBack() {
  theoretical_peak_set_.Clear();
  Peptide* peptide = queue_.back();
  peptide->ComputeTheoreticalPeaks(&theoretical_peak_set_, dia_mode_);
  ion_inverted_index_.insert_peaks(peptide);
}

int ActivePeptideQueue::SetActiveRange(double min_range, double max_range, double max_exp_peak_mz) {

  // Reset the Tailor score histogram for the new spectrum
  ResetTailorHistogram();

  // min_range and max_range have been introduced to fix a bug
  // introduced by m/z selection. see #222 in sourceforge
  // this has to be true:
  // min_range <= min_mass <= max_mass <= max_range

  // queue front() is lightest; back() is heaviest

  // delete anything already loaded that falls below min_range
  while (!queue_.empty() && queue_.front()->Mass() < min_range) {
    Peptide* peptide = queue_.front();
    queue_.pop_front();
    ion_inverted_index_.pop_peaks(peptide);
    delete peptide;
  }
  total_peptides_loaded_ = 0;
  candidate_peptide_count_ = 0;
  candidate_target_count_ = 0;
  candidate_decoy_count_ = 0;
  theoretical_peak_set_.max_exp_peak_mz_ = max_exp_peak_mz + 100.00;

  // Enqueue all peptides that are not yet queued but are lighter than
  // max_range. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.
  bool done = false;
  //Modified for tailor score calibration method by AKF
  if (queue_.empty() || queue_.back()->Mass() <= max_range || queue_.size() < min_candidates_) {
    if (!queue_.empty()) {
      ComputeTheoreticalPeaksBack();
    }
    while (!(done = reader_->Done())) {
      // read all peptides lighter than max_range
      reader_->Read(&current_pb_peptide_);
      if (current_pb_peptide_.mass() < min_range) {
        // we would delete current_pb_peptide_;
        continue; // skip peptides that fall below min_range
      }
      Peptide* peptide = new Peptide(current_pb_peptide_, proteins_, curScoreFunction_, locations_);
      assert(peptide != NULL);
      queue_.push_back(peptide);
      // ion_inv_.insert(peptide);
      //Modified for tailor score calibration method by AKF
      if (peptide->Mass() > max_range && queue_.size() > min_candidates_) {
        break;
      }
      ComputeTheoreticalPeaksBack();
    }
  }
  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!queue_.empty() || done);

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  if (queue_.empty()) {
    return 0;
  }

  total_peptides_loaded_ = 0;
  begin_ = queue_.begin();
  end_ = queue_.end();
  return queue_.size();
}

// Reset state: clear E-value histogram and regression parameters before scoring a spectrum
void ActivePeptideQueue::ResetEValueHistogram() {
  memset(evalue_histogram_, 0, sizeof(evalue_histogram_));
  evalue_score_count_ = 0;
  evalue_slope_ = 0.0;
  evalue_intercept_ = 0.0;
  evalue_reg_start_bin_ = 0;
  evalue_reg_end_bin_ = 0;
  evalue_max_bin_ = 0;
}

// Add a score to the E-value histogram (collects ALL scores, target + decoy)
void ActivePeptideQueue::AddToEValueHistogram(double score) {
  int bin = ScoreToBin(score);
  evalue_histogram_[bin]++;
  evalue_score_count_++;
}

// Map a score to an E-value histogram bin (score × 10, clamped to [0, EVALUE_HISTOGRAM_SIZE))
int ActivePeptideQueue::ScoreToBin(const double score) const {
  int bin = static_cast<int>(score * 10.0 + 0.5);
  if (bin < 0) bin = 0;
  if (bin >= EVALUE_HISTOGRAM_SIZE) bin = EVALUE_HISTOGRAM_SIZE - 1;
  return bin;
}

// Fit regression model on the E-value score histogram; enables ComputeEValue for this spectrum.
// For XCorr: fits log10(cumulative count) vs bin index.
// For HyperScore: fits log10(survival function S(x)) vs bin index.
void ActivePeptideQueue::FitEValueRegression() {
  if (evalue_score_count_ < MIN_SCORES_FOR_REGRESSION) {
    evalue_slope_ = 0.0;
    evalue_intercept_ = 0.0;
    return;
  }

  if (curScoreFunction_ == HYPERSCORE) {
    double rsq;
    LinearRegressionHyperScore(evalue_histogram_, EVALUE_HISTOGRAM_SIZE, &evalue_slope_, &evalue_intercept_, &rsq);
  } else {
    LinearRegression(evalue_histogram_, &evalue_slope_, &evalue_intercept_,
                     &evalue_max_bin_, &evalue_reg_start_bin_, &evalue_reg_end_bin_);
  }
}

double ActivePeptideQueue::ComputeEValue(double score) const {
  if (evalue_score_count_ < MIN_SCORES_FOR_REGRESSION) return EVALUE_UPPER_BOUND;
  if (evalue_slope_ == 0.0 && evalue_intercept_ == 0.0) return EVALUE_UPPER_BOUND;

  int bin = ScoreToBin(score);
  double exponent = evalue_slope_ * bin + evalue_intercept_;
  double eval = pow(10.0, exponent);
  if (eval > EVALUE_UPPER_BOUND) eval = EVALUE_UPPER_BOUND;
  if (eval < 0.0) eval = 0.0;
  return eval;
}

void ActivePeptideQueue::LinearRegression(const int* histogram, double* slope, double* intercept,
                             int* max_bin, int* start_bin, int* end_bin) {
  int i;
  for (i = EVALUE_HISTOGRAM_SIZE - 1; i >= 0; --i) {
    if (histogram[i] > 0) break;
  }
  *max_bin = i;

  if (*max_bin < 0) {
    *slope = 0.0;
    *intercept = 0.0;
    *start_bin = -1;
    *end_bin = -1;
    return;
  }

  double totalScores = 0.0;
  for (i = 0; i <= *max_bin; ++i) totalScores += histogram[i];
  if (totalScores < MIN_SCORES_FOR_REGRESSION) {
    *slope = 0.0;
    *intercept = 0.0;
    *start_bin = -1;
    *end_bin = -1;
    return;
  }

  int iEndBin = 0;
  bool foundFirstNonZero = false;

  for (i = 0; i < *max_bin; ++i) {
    if (histogram[i] == 0 && foundFirstNonZero && i >= 10) {
      if (histogram[i + 1] == 0 || i + 1 == *max_bin) {
        iEndBin = (i > 0) ? i - 1 : i;
        break;
      }
    }
    if (histogram[i] != 0) foundFirstNonZero = true;
  }

  if (i == *max_bin) {
    iEndBin = *max_bin;
    if (*max_bin >= 10) {
      for (i = *max_bin; i >= *max_bin - 5; --i) {
        if (histogram[i] == 0) {
          iEndBin = i;
          if (*max_bin <= 20) break;
        }
      }
      if (iEndBin == *max_bin) iEndBin = *max_bin - 1;
    }
  }
  *end_bin = iEndBin;

  double cumulative[EVALUE_HISTOGRAM_SIZE];
  cumulative[*end_bin] = histogram[*end_bin];
  for (i = *end_bin - 1; i >= 0; --i) {
    cumulative[i] = cumulative[i + 1] + histogram[i];
  }

  double logCumulative[EVALUE_HISTOGRAM_SIZE];
  for (i = 0; i <= *end_bin; ++i) {
    if (cumulative[i] > 0.0)
      logCumulative[i] = log10(cumulative[i]);
    else
      logCumulative[i] = 0.0;
  }

  int iStart = *end_bin - 5;
  int numZeros = 0;
  for (i = iStart; i <= *end_bin; ++i)
    if (logCumulative[i] == 0.0) numZeros++;
  iStart -= numZeros;
  if (iStart < 0) iStart = 0;

  double Mx, My, a, b;
  Mx = My = a = b = 0.0;

  while (iStart >= 0 && *end_bin > iStart + 2) {
    double Sx = 0.0, Sxy = 0.0, SumX = 0.0, SumY = 0.0;
    int n = 0;

    for (i = iStart; i <= *end_bin; ++i) {
      if (cumulative[i] > 0.0) {
        SumX += i;
        SumY += logCumulative[i];
        ++n;
      }
    }

    if (n > 0) {
      Mx = SumX / n;
      My = SumY / n;
    } else {
      Mx = My = 0.0;
    }

    for (i = iStart; i <= *end_bin; ++i) {
      if (cumulative[i] > 0.0) {
        double dX = i - Mx;
        double dY = logCumulative[i] - My;
        Sx += dX * dX;
        Sxy += dX * dY;
      }
    }

    if (Sx > 0.0)
      b = Sxy / Sx;
    else
      b = 0.0;

    if (b < 0.0)
      break;
    else
      iStart--;
  }

  a = My - b * Mx;
  *start_bin = iStart;
  *slope = b;
  *intercept = a;
}

// hyperscore_expect.cpp
//
// Reads hyperscore histogram .txt files and computes an E-value for each
// spectrum by fitting a linear regression to the log10-transformed
// survival function in the histogram tail.
//
// File format (per input file):
//   Line 1 : <ref_slope>\t<ref_intercept>   (pre-computed reference params; skipped)
//   Lines 2+: <score>\t<count>              (integer hyperscore 0..1999)
//
// Algorithm:
//   1. Build S(x) = P(score >= x) = cumulative_count_from_right / total
//   2. Identify the tail window using two complementary guards:
//      a. Gap detection: stop at the first run of GAP_THRESHOLD consecutive
//         zero bins, suppressing isolated high-score outlier clusters.
//      b. Survival thresholds: only fit scores where
//            S_min <= S(x) <= S_max
//         S_max = SURV_MAX (default 0.10) -- skip the bulk region (S > 10%)
//                where the distribution is not yet in exponential decay.
//         S_min = max(5/N, 1e-3)          -- skip the sparse floor where
//                individual counts are too few for reliable log estimates.
//         If fewer than MIN_TAIL_POINTS survive the S_max cut, the filter
//         is relaxed (S_max removed) and only S_min is applied.
//   3. Fit   log10(S(x)) = slope * x + intercept   by ordinary least squares.
//   4. E-value at score s:  E = 10^(slope * s + intercept)
//
// Output (stdout, tab-aligned):
//   filename  evalue  slope  intercept  r2  top_score  surv_count  ref_evalue  ref_slope  ref_intercept
//
// E-value units: expected number of random candidates scoring >= top_score.
//   evalue     = N * 10^(slope * top_score + intercept)
//                where N = total histogram count and intercept is from the
//                fraction fit, so the result is an absolute count (>1 = poor match).
//   surv_count = raw histogram count at scores >= top_score (ground truth).
//   ref_evalue = 10^(ref_slope * top_score + ref_intercept), the pre-computed
//                reference E-value already in absolute-count units.

bool ActivePeptideQueue::LinearRegressionHyperScore(const int* histogram,
                     int hist_size,
                     double* slope,
                     double* intercept,
                     double* rsq)
{
   // --- total counts ---
   long long llTotal = 0;
   for (int i = 0; i < hist_size; i++)
      llTotal += histogram[i];
   if (llTotal == 0)
      return false;

   // --- mode: score with the highest count ---
   int iMode = 0;
   for (int i = 1; i < hist_size; i++)
   {
      if (histogram[i] > histogram[iMode])
         iMode = i;
   }

   // --- survival fraction S(x) = count(score >= x) / total ---
   std::vector<double> vdS(hist_size, 0.0);
   {
      long long llCum = 0;
      for (int i = hist_size - 1; i >= 0; i--)
      {
         llCum   += histogram[i];
         vdS[i]   = (double)llCum / (double)llTotal;
      }
   }

   // --- tail region: [iMode+1, iTailEnd] ---
   // Walk rightward from the mode; stop at the first run of GAP_THRESHOLD
   // consecutive zero bins.  This excludes isolated high-score outliers
   // (e.g., a handful of true identifications far above the noise bulk)
   // that would otherwise flatten the regression slope.
   int iTailEnd  = iMode;
   int iGapCount = 0;

   for (int i = iMode + 1; i < hist_size; i++)
   {
      if (histogram[i] > 0)
      {
         iTailEnd  = i;
         iGapCount = 0;
      }
      else
      {
         ++iGapCount;
         if (iGapCount >= GAP_THRESHOLD)
            break;
      }
   }

   int iTailStart = iMode + 1;
   if (iTailEnd < iTailStart)
      return false;   // nothing above the mode within the gap threshold

   // --- survival bounds ---
   // S_min: at least 5 counts expected, or 0.1% of N (whichever is larger).
   // S_max: SURV_MAX -- skip the bulk region where S(x) is still high.
   const double dSurvMin = fmax(5.0 / (double)llTotal, 1e-3);
   const double dSurvMax = SURV_MAX;

   // --- collect (x, log10(S(x))) pairs: primary pass with both bounds ---
   std::vector<double> vdX;
   std::vector<double> vdY;
   vdX.reserve(iTailEnd - iTailStart + 1);
   vdY.reserve(iTailEnd - iTailStart + 1);

   for (int i = iTailStart; i <= iTailEnd; i++)
   {
      if (vdS[i] >= dSurvMin && vdS[i] <= dSurvMax)
      {
         vdX.push_back((double)i);
         vdY.push_back(log10(vdS[i]));
      }
   }

   // --- fallback: if S_max filter leaves too few points, drop it ---
   // This handles distributions where S never falls below SURV_MAX in the
   // gap-truncated window (e.g., very sparse histograms or flat tails).
   if ((int)vdX.size() < MIN_TAIL_POINTS)
   {
      vdX.clear();
      vdY.clear();
      for (int i = iTailStart; i <= iTailEnd; i++)
      {
         if (vdS[i] >= dSurvMin)
         {
            vdX.push_back((double)i);
            vdY.push_back(log10(vdS[i]));
         }
      }
   }

   int iN = (int)vdX.size();
   if (iN < MIN_TAIL_POINTS)
      return false;

   // --- ordinary least squares: y = slope * x + intercept ---
   double dSumX  = 0.0;
   double dSumY  = 0.0;
   double dSumXX = 0.0;
   double dSumXY = 0.0;

   for (int i = 0; i < iN; i++)
   {
      dSumX  += vdX[i];
      dSumY  += vdY[i];
      dSumXX += vdX[i] * vdX[i];
      dSumXY += vdX[i] * vdY[i];
   }

   double dDenom = (double)iN * dSumXX - dSumX * dSumX;
   if (fabs(dDenom) < 1e-12)
      return false;

   *slope     = ((double)iN * dSumXY - dSumX * dSumY) / dDenom;
   *intercept = (dSumY - *slope * dSumX) / (double)iN;

   // --- R^2: fraction of log10(S) variance explained by the linear fit ---
   double dMeanY    = dSumY / (double)iN;
   double dSSTot    = 0.0;
   double dSSRes    = 0.0;
   for (int i = 0; i < iN; i++)
   {
      double dYhat  = *slope * vdX[i] + *intercept;
      double dDev   = vdY[i] - dMeanY;
      double dRes   = vdY[i] - dYhat;
      dSSTot += dDev * dDev;
      dSSRes += dRes * dRes;
   }
   *rsq = (dSSTot < 1e-15) ? 1.0 : 1.0 - dSSRes / dSSTot;

   return true;
}
