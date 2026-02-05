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

}



