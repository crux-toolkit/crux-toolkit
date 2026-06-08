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
  min_candidates_ = 3000;
  nPeptides_ = 0;
  nCandPeptides_ = 0;
  CandPeptidesTarget_ = 0;
  CandPeptidesDecoy_ = 0;  


  score_histogram_offset_= 10;   // This should be ok for XCorr. It does not matter for HyperScore; TAILOR_OFFSET
  score_scale_factor_ = 100;   // TODO: may be fine-tuned experimentally.
  max_score_ = 100;
  highest_bin_ = 0;
  score_count_ = 1; // so that the score histogram vector will be set to zero.
  score_histogram_.resize((max_score_+score_histogram_offset_)*score_scale_factor_); // The maximum score can be 100. Can be increased
  ResetHist();
  curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function")); 

  // e-value fields initialization
  memset(xcorrHistogram_, 0, sizeof(xcorrHistogram_));
  decoyCount_ = 0;
  slope_ = 0.0;
  intercept_ = 0.0;
  startCorr_ = 0;
  nextCorr_ = 0;
  maxCorr_ = 0;
}

void ActivePeptideQueue::AddScoreToHist(double score, int match_cnt) {
  int bin = round( (score + score_histogram_offset_)*score_scale_factor_ );
  score_histogram_[bin] += 1;
  ++score_count_;
  total_match_ += match_cnt;
  if (highest_bin_ < bin)
    highest_bin_ = bin;
}

void ActivePeptideQueue::ResetHist() {
  total_match_ = 0;
  highest_bin_ = 0;
  if (score_count_ == 0)
    return;
  //memset(score_histogram_.data(), 0, score_histogram_.capacity() * sizeof(int));
  std::fill(score_histogram_.begin(), score_histogram_.end(), 0); // reset the values of the score histogram    
  score_count_ = 0;
}

double ActivePeptideQueue::getTailorQuantileFromHistogram(){
  int bin = highest_bin_+1; //score_histogram_.capacity()-1;
  int count = 0;
  double quantile_score_ = 1.0;

  int quantile_pos = (int)(TAILOR_QUANTILE_TH*(double)score_count_+0.5)-1; // zero indexed
  if (quantile_pos < 2) 
    quantile_pos = 2;  // the third element
  if (quantile_pos >= score_count_) 
    quantile_pos = score_count_-1; // the last element

  while (count < quantile_pos ) {
    count += score_histogram_[bin--];
  }
  ++bin;
  return (double)bin/score_scale_factor_ - score_histogram_offset_; // Make sure scores positive in case of XCorr;
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

  // Reset the score histograms
  ResetHist();

  //min_range and max_range have been introduced to fix a bug
  //introduced by m/z selection. see #222 in sourceforge
  //this has to be true:
  // min_range <= min_mass <= max_mass <= max_range

  // queue front() is lightest; back() is heaviest

  // delete anything already loaded that falls below min_range
  while (!queue_.empty() && queue_.front()->Mass() < min_range) {
    Peptide* peptide = queue_.front();
    queue_.pop_front();
    ion_inverted_index_.pop_peaks(peptide);
    delete peptide;
  }
  nPeptides_ = 0;
  nCandPeptides_ = 0;
  CandPeptidesTarget_ = 0;
  CandPeptidesDecoy_ = 0;
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

  nPeptides_ = 0;
  begin_ = queue_.begin();  
  end_ = queue_.end();
  return queue_.size();

}// reset state for a new spectrum
void ActivePeptideQueue::BeginSpectrum() {
  memset(xcorrHistogram_, 0, sizeof(xcorrHistogram_));
  decoyCount_ = 0;
  slope_ = 0.0;
  intercept_ = 0.0;
  startCorr_ = 0;
  nextCorr_ = 0;
  maxCorr_ = 0;
}

// add XCorr of a decoy match to the histogram of the current spectrum
void ActivePeptideQueue::AddDecoyXCorr(double xcorr) {
  // int bin = static_cast<int>(xcorr * 10.0 + 0.5);
  // if (bin < 0) bin = 0;
  // if (bin >= HISTO_SIZE) bin = HISTO_SIZE - 1;
  int bin = XCorrToBin(xcorr);   
  xcorrHistogram_[bin]++;
  decoyCount_++;
}

void ActivePeptideQueue::EndSpectrum() {
  if (decoyCount_ < MIN_DECOY_COUNT) {
    slope_ = 0.0;
    intercept_ = 0.0;
    return;
  }

  if (curScoreFunction_ == HYPERSCORE) {
    double pdRsq;
    LinearRegressionHyperScore(xcorrHistogram_, HISTO_SIZE, &slope_, &intercept_, &pdRsq);
  } else {
    LinearRegression(xcorrHistogram_, &slope_, &intercept_, &maxCorr_,
                    &startCorr_, &nextCorr_);
  }

  // LinearRegressionHyperScore(const int* piHistogram, int iHistSize, double* pdSlope, double* pdIntercept, double* pdRsq);
                   
  // FILE * fp = fopen("score_histogram_scanID.txt", "w");
  // fprintf(fp, "%lf\t%lf\n" ,slope_, intercept_);
  // for (int i = 0; i< HISTO_SIZE; ++i) {
  //   fprintf(fp, "%d\t%d\n", i,xcorrHistogram_[i]);
  // }
  // fclose(fp);
}

double ActivePeptideQueue::ComputeEValue(double xcorr) const {
  if (decoyCount_ < MIN_DECOY_COUNT) return MAX_EVALUE;
  if (slope_ == 0.0 && intercept_ == 0.0) return MAX_EVALUE;

  int bin = XCorrToBin(xcorr); 
  double exponent = slope_ * bin + intercept_;
  double eval = pow(10.0, exponent);
  if (eval > MAX_EVALUE) eval = MAX_EVALUE;
  if (eval < 0.0) eval = 0.0;
  return eval;
}

int ActivePeptideQueue::XCorrToBin(const double xcorr) const {
  int bin = static_cast<int>(xcorr * 10.0 + 0.5);
  if (bin < 0) bin = 0;
  if (bin >= HISTO_SIZE) bin = HISTO_SIZE - 1;
  return bin;
}

void ActivePeptideQueue::LinearRegression(int* histogram, double* slope, double* intercept,
                             int* maxCorr, int* startCorr, int* nextCorr) { 
  int i;
  for (i = HISTO_SIZE - 1; i >= 0; --i) {
    if (histogram[i] > 0) break;
  }
  *maxCorr = i;

  if (*maxCorr < 0) {
    *slope = 0.0;
    *intercept = 0.0;
    *startCorr = -1;
    *nextCorr = -1;
    return;
  }

  double totalDecoys = 0.0;
  for (i = 0; i <= *maxCorr; ++i) totalDecoys += histogram[i];
  if (totalDecoys < MIN_DECOY_COUNT) {
    *slope = 0.0;
    *intercept = 0.0;
    *startCorr = -1;
    *nextCorr = -1;
    return;
  }

  int iNextCorr = 0;
  bool foundFirstNonZero = false;

  for (i = 0; i < *maxCorr; ++i) {
    if (histogram[i] == 0 && foundFirstNonZero && i >= 10) {
      if (histogram[i + 1] == 0 || i + 1 == *maxCorr) {
        iNextCorr = (i > 0) ? i - 1 : i;
        break;
      }
    }
    if (histogram[i] != 0) foundFirstNonZero = true;
  }

  if (i == *maxCorr) {
    iNextCorr = *maxCorr;
    if (*maxCorr >= 10) {
      for (i = *maxCorr; i >= *maxCorr - 5; --i) {
        if (histogram[i] == 0) {
          iNextCorr = i;
          if (*maxCorr <= 20) break;
        }
      }
      if (iNextCorr == *maxCorr) iNextCorr = *maxCorr - 1;
    }
  }
  *nextCorr = iNextCorr;

  double cumulative[HISTO_SIZE];
  cumulative[*nextCorr] = histogram[*nextCorr];
  for (i = *nextCorr - 1; i >= 0; --i) {
    cumulative[i] = cumulative[i + 1] + histogram[i];
  }

  double logCumulative[HISTO_SIZE];
  for (i = 0; i <= *nextCorr; ++i) {
    if (cumulative[i] > 0.0)
      logCumulative[i] = log10(cumulative[i]);
    else
      logCumulative[i] = 0.0;
  }

  int iStart = *nextCorr - 5;
  int numZeros = 0;
  for (i = iStart; i <= *nextCorr; ++i)
    if (logCumulative[i] == 0.0) numZeros++;
  iStart -= numZeros;
  if (iStart < 0) iStart = 0;

  double Mx, My, a, b;
  Mx = My = a = b = 0.0;

  while (iStart >= 0 && *nextCorr > iStart + 2) {
    double Sx = 0.0, Sxy = 0.0, SumX = 0.0, SumY = 0.0;
    int n = 0;

    for (i = iStart; i <= *nextCorr; ++i) {
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

    for (i = iStart; i <= *nextCorr; ++i) {
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
  *startCorr = iStart;
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

bool ActivePeptideQueue::LinearRegressionHyperScore(const int* piHistogram,
                     int        iHistSize,
                     double*    pdSlope,
                     double*    pdIntercept,
                     double*    pdRsq)
{
   // --- total counts ---
   long long llTotal = 0;
   for (int i = 0; i < iHistSize; i++)
      llTotal += piHistogram[i];
   if (llTotal == 0)
      return false;

   // --- mode: score with the highest count ---
   int iMode = 0;
   for (int i = 1; i < iHistSize; i++)
   {
      if (piHistogram[i] > piHistogram[iMode])
         iMode = i;
   }

   // --- survival fraction S(x) = count(score >= x) / total ---
   std::vector<double> vdS(iHistSize, 0.0);
   {
      long long llCum = 0;
      for (int i = iHistSize - 1; i >= 0; i--)
      {
         llCum   += piHistogram[i];
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

   for (int i = iMode + 1; i < iHistSize; i++)
   {
      if (piHistogram[i] > 0)
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

   *pdSlope     = ((double)iN * dSumXY - dSumX * dSumY) / dDenom;
   *pdIntercept = (dSumY - *pdSlope * dSumX) / (double)iN;

   // --- R^2: fraction of log10(S) variance explained by the linear fit ---
   double dMeanY    = dSumY / (double)iN;
   double dSSTot    = 0.0;
   double dSSRes    = 0.0;
   for (int i = 0; i < iN; i++)
   {
      double dYhat  = *pdSlope * vdX[i] + *pdIntercept;
      double dDev   = vdY[i] - dMeanY;
      double dRes   = vdY[i] - dYhat;
      dSSTot += dDev * dDev;
      dSSRes += dRes * dRes;
   }
   *pdRsq = (dSSTot < 1e-15) ? 1.0 : 1.0 - dSSRes / dSSTot;

   return true;
}
