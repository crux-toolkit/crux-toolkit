// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
#include <deque>
#include <gflags/gflags.h>
#include "records.h"
#include "peptides.pb.h"
#include "peptide.h"
#include "PeptideDiskLoader.h"
#include "records_to_vector-inl.h"
#include "theoretical_peak_set.h"
#include "compiler.h"
#include "app/TideMatchSet.h"
#include <map> 
#include <limits>
#include <iostream>
#define CHECK(x) GOOGLE_CHECK((x))

PeptideDiskLoader::PeptideDiskLoader(RecordReader* reader,
                                       const vector<const pb::Protein*>& proteins, 
                                       vector<const pb::AuxLocation*>* locations, 
                                       bool dia_mode, int thread_num)
  : reader_(reader),
    proteins_(proteins),
    locations_(locations),
    dia_mode_(dia_mode),
    windows_(thread_num, nullptr) {

  CHECK(reader_->OK());
  min_candidates_ = 30;
  for (size_t i = 0; i < windows_.size(); ++i) {
    windows_[i] = new RollingPeptideWindow(this);
  }
}

PeptideDiskLoader::~PeptideDiskLoader() {
  for (size_t i = 0; i < windows_.size(); ++i) {
    delete windows_[i];
  }
}
/*
void PeptideDiskLoader::ComputeTheoreticalPeak(size_t i) {
  static thread_local TheoreticalPeakSetBYSparse theoretical_peak_set(1000); // Probably overkill, but no harm
  theoretical_peak_set.Clear();
  Peptide* peptide = getPeptide(i);
  peptide->Activate(&theoretical_peak_set, dia_mode_);
}

void PeptideDiskLoader::ComputeTheoreticalPeak(Peptide* peptide) {
  static thread_local TheoreticalPeakSetBYSparse theoretical_peak_set(1000); // Probably overkill, but no harm
  theoretical_peak_set.Clear();
  peptide->Activate(&theoretical_peak_set, dia_mode_);
}
*/


bool PeptideDiskLoader::pushBack(RollingPeptideWindow* window) {
  if (window->last_ != nullptr && window->last_->next.load() != nullptr) {
    window->window_.push_back(window->last_->next.load());
    window->end_++;
    return true;
  }
  std::lock_guard<SpinLock> lock(m_);
  assert(window->end_ <= end_);
  if (window->end_ < end_) {
    window->end_++;
    window->window_.push_back(getPeptideUnsafe(window->end_ - 1));
    return true;
  }
  if (!reader_->Done()) {
    reader_->Read(&current_pb_peptide_);
    Peptide* peptide = new Peptide(current_pb_peptide_, proteins_, locations_);
    assert(peptide != NULL);
    queue_.push_back(peptide);
    end_++;
    window->end_++;
    // ComputeTheoreticalPeak(peptide);
    window->window_.push_back(getPeptideUnsafe(window->end_ - 1));
    return true;
  }
  return false;
}

bool PeptideDiskLoader::popFront(RollingPeptideWindow* window) {
  if (window->begin_ == window->end_) {
    return false;
  }
  assert(window->window_.size() > 0);
  Peptide* peptide = window->window_.front();
  assert(peptide);
  if (peptide->Drop() == windows_.size() - 1) {
    std::lock_guard<SpinLock> lock(m_);
    queue_.pop_front();
    delete peptide;
    begin_++;
    window->last_ = nullptr;
  }
  window->window_.pop_front();
  window->begin_++;
  return true;
}

Peptide* PeptideDiskLoader::getPeptideUnsafe(size_t i) {
  Peptide* peptide = queue_.at(i - begin_);
  return peptide;
}
/*
Peptide* PeptideDiskLoader::getComputedPeptide(size_t i) {
  ComputeTheoreticalPeak(i);
  return getPeptide(i);
}
*/

const std::vector<RollingPeptideWindow*> PeptideDiskLoader::GetRollingPeptideWindows() const {
  return windows_;
}

RollingPeptideWindow::RollingPeptideWindow(PeptideDiskLoader* queue) : 
  theoretical_peak_set_(1000),
  queue_(queue) {
}

void RollingPeptideWindow::ActivatePeptides(bool dia_mode) {
  for (size_t i = begin(); i < end(); ++i) {
    Peptide* peptide = GetPeptide(i);
    peptide->TryActivate(&theoretical_peak_set_, dia_mode);
  }
}

bool RollingPeptideWindow::PushBack() {
  if (!queue_->pushBack(this)) {
    return false;
  }
  assert(!window_.empty());
  last_ = window_.back();
  window_.pop_back();
  if (!window_.empty()) {
    window_.back()->next.store(last_);
  }
  window_.push_back(last_);
  return true;
}

bool RollingPeptideWindow::PopFront() {
  bool moved = queue_->popFront(this);
  return moved;
}

int RollingPeptideWindow::SetActiveRange(double min_range, double max_range, vector<double>* min_mass, vector<double>* max_mass) {
  // queue front() is lightest; back() is heaviest

  nCandPeptides_ = 0;
  CandPeptidesTarget_ = 0;
  CandPeptidesDecoy_ = 0;

  // delete anything already loaded that falls below min_range
  bool done = false;
  while (this->size() < queue_->min_candidates_ || this->front()->Mass() < min_range || this->back()->Mass() <= max_range) {
    if (!this->empty() && this->front()->Mass() < min_range) { // should move front end
      this->PopFront();
      continue;
    }
    done = this->PushBack();
    if (!done) {
      break;
    }

  }

  /*
  while (!queue_.empty() && queue_.front()->Mass() < min_range) {
    Peptide* peptide = queue_.front();
    queue_.pop_front();
    delete peptide;
  }
  */
  nPeptides_ = size();

  // Enqueue all peptides that are not yet queued but are lighter than
  // max_range. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.

  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  
  // assert(!empty() || done);  AKF commented out for now. REthing if we need this.

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  if (empty()) {
    return 0;
  }
  min_mass_ = min_mass;
  max_mass_ = max_mass;

  return nPeptides_;
}