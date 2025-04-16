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
    windows_(thread_num) {

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
// bool PeptideDiskLoader::isWithinIsotope(vector<double>* min_mass, vector<double>* max_mass, double mass, int* isotope_idx) {
//   for (int i = *isotope_idx; i < min_mass->size(); ++i) {
//     if (mass >= (*min_mass)[i] && mass <= (*max_mass)[i]) {
//       if (i > *isotope_idx) {
//         *isotope_idx = i;
//       }
//       return true;
//     }
//   }
//   return false;
// }

bool PeptideDiskLoader::pushBack(RollingPeptideWindow* window) {
  std::lock_guard<boost::shared_mutex> lock(m_);
  assert(window->end_ <= end_);
  if (window->end_ < end_) {
    window->end_++;
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
    return true;
  }
  return false;
}

bool PeptideDiskLoader::popFront(RollingPeptideWindow* window) {
  std::lock_guard<boost::shared_mutex> lock(m_);
  if (window->begin_ == window->end_) {
    return false;
  }
  Peptide* peptide = getPeptideUnsafe(window->begin_);
  peptide->Drop();

  if (peptide->GetDropsCounter() == windows_.size()) {
    delete peptide;
    queue_.pop_front();
    begin_++;
  }
  window->begin_++;
  return true;
}

Peptide* PeptideDiskLoader::getPeptide(size_t i) {
  m_.lock_shared();
  Peptide* peptide = queue_.at(i - begin_);
  m_.unlock_shared();
  return peptide;
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

bool RollingPeptideWindow::PushBack() {
  return queue_->pushBack(this);
}

bool RollingPeptideWindow::PopFront() {
  bool moved = queue_->popFront(this);
  return moved;
}

int RollingPeptideWindow::SetActiveRange(double min_range, double max_range) {
  // queue front() is lightest; back() is heaviest

  // delete anything already loaded that falls below min_range
  // std::cout << "\n";
  bool done = false;
  while (this->size() < queue_->min_candidates_ || this->front()->Mass() < min_range || this->back()->Mass() <= max_range) {
    if (!this->empty() && this->front()->Mass() < min_range) { // should move front end
      this->PopFront();
      continue;
    }
    if (!(done = this->PushBack())) {
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
  nCandPeptides_ = size();
  CandPeptidesTarget_ = 0;
  CandPeptidesDecoy_ = 0;

  // Enqueue all peptides that are not yet queued but are lighter than
  // max_range. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.

  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!empty() || done);

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  if (empty()) {
    return 0;
  }

  return nCandPeptides_;
}