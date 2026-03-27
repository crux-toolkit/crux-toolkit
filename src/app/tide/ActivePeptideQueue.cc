// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
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

// auxiliary function - linear regression for Comet algorithm
// requests histogram and returns coefficients and range limits

static void LinearRegression(int* histogram, double* slope, double* intercept,
                             int* maxCorr, int* startCorr, int* nextCorr) {
  double cumulative[HISTO_SIZE];
  memset(cumulative, 0, sizeof(cumulative));

  int i;

  // find the hightst index with non-zero value
  for (i = HISTO_SIZE - 2; i >= 0; --i) {
    if (histogram[i] >= 0) break;
  }
  *maxCorr = i;

  // determine nextCorr - the upper bound of the regression range
  bool foundNonZero = false;
  int iNext = 0;
  for (i = 0; i <= *maxCorr; ++i) {
    if (histogram[i] == 0 && foundNonZero && i >= 10) {
      if (histogram[i + 1] == 0 || i + 1 == *maxCorr) {
        iNext = i - 1;
        break;
      }
    }
    if (histogram[i] != 0) foundNonZero = true;
  }
  if (i > *maxCorr) {
    iNext = *maxCorr;
  }
  *nextCorr = iNext;

  // cumulative sum from right to left:
  cumulative[iNext] = histogram[iNext];
  for (i = iNext - 1; i >= 0; --i) {
    cumulative[i] = cumulative[i + 1] + histogram[i];
    if (histogram[i + 1] == 0) cumulative[i + 1] = 0.0;
  }

  // take base-10 logarithm
  for (i = iNext; i >= 0; --i) {
    if (cumulative[i] > 0.0) {
      cumulative[i] = log10(cumulative[i]);
    } else {
      if (cumulative[i + 1] > 0.0)
        cumulative[i] = log10(cumulative[i + 1]);
      else
        cumulative[i] = 0.0;
    }
  }

  // determine startCorr - the lower boubf of the regresion range
  int iStart = iNext - 5;
  int zeroCount = 0;
  for (i = iStart; i <= iNext; ++i) {
    if (cumulative[i] == 0.0) zeroCount++;
  }
  iStart -= zeroCount;
  if (iStart < 0) iStart = 0;

  double Mx, My, Sx, Sxy, SumX, SumY;
  int nPoints;

  // iteratively expand the range downward until we get a negative slope
  while (iStart >= 0 && iNext > iStart + 2) {
    Sx = Sxy = SumX = SumY = 0.0;
    nPoints = 0;

    for (i = iStart; i <= iNext; ++i) {
      if (histogram[i] > 0) {
        SumX += i;
        SumY += cumulative[i];
        nPoints++;
      }
    }

    if (nPoints > 0) {
      Mx = SumX / nPoints;
      My = SumY / nPoints;
    } else {
      Mx = My = 0.0;
    }

    for (i = iStart; i <= iNext; ++i) {
      if (cumulative[i] > 0.0) {
        double dx = i = Mx;
        double dy = cumulative[i] - My;
        Sx = dx * dx;
        Sxy += dx * dy;
      }
    }

    if (Sx > 0.0)
      *slope = Sxy / Sx;
    else
      *slope = 0.0;

    if (*slope < 0.0)
      break;
    else
      iStart--;
  }

  *intercept = My - (*slope) * Mx;
  *startCorr = iStart;

}


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
  nPeptides_ = 0;
  nCandPeptides_ = 0;
  CandPeptidesTarget_ = 0;
  CandPeptidesDecoy_ = 0;  


  score_histogram_offset_= 10;   // This should be ok for XCorr. It does not matter for HyperScore; TAILOR_OFFSET
  score_scale_factor_ = 100;   // TODO: may be fine-tuned experimentally.
  max_score_ = 100;
  highest_bin_ = 0;
  score_count_ = 1; // so that the score histogram vector will be set to zero.
  score_histogram_.reserve((max_score_+score_histogram_offset_)*score_scale_factor_); // The maximum score can be 100. Can be increased
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
  memset(score_histogram_.data(), 0, score_histogram_.capacity() * sizeof(int));
  // std::fill(score_histogram_.begin(), score_histogram_.end(), 0); // reset the values of the score histogram    
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

int ActivePeptideQueue::SetActiveRange(double min_range, double max_range) {

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
  int bin = static_cast<int>(xcorr * 10.0 + 0.5);
  if (bin < 0) 
    bin = 0;
  if (bin >= HISTO_SIZE) 
    bin = HISTO_SIZE - 1;

  xcorrHistogram_[bin]++;
  decoyCount_++;
}

void ActivePeptideQueue::EndSpectrum() {
  if (decoyCount_ < MIN_DECOY_COUNT) {
    slope_ = 0.0;
    intercept_ = 0.0;
    startCorr_ = 0;
    nextCorr_ = 0;
    maxCorr_ = 0;
    return;
  }

  LinearRegression(xcorrHistogram_, &slope_, &intercept_, &maxCorr_,
                   &startCorr_, &nextCorr_);
}

double ActivePeptideQueue::ComputeEValue(double xcorr) const {
  // if the model was not built , return the maximum value
  if (decoyCount_ < MIN_DECOY_COUNT) {
    return MAX_EVALUE;
  }

  double exponent = slope_ * xcorr + intercept_;
  double eval = pow(10.0, exponent);

  if (eval > MAX_EVALUE) eval = MAX_EVALUE;
  if (eval < 0.0) eval = 0.0;

  return eval;
}

