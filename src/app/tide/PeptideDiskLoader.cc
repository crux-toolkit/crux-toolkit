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
#define CHECK(x) GOOGLE_CHECK((x))

RollingPeptideWindow::RollingPeptideWindow(PeptideDiskLoader* queue) : queue_(queue) {
}

bool RollingPeptideWindow::PushBack() {
  return queue_->pushBack(this);
}

bool RollingPeptideWindow::PopFront() {
  bool moved = queue_->popFront(this);
  return moved;
}

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
}

void PeptideDiskLoader::ComputeTheoreticalPeak(size_t i) {
  static thread_local TheoreticalPeakSetBYSparse theoretical_peak_set(1000); // Probably overkill, but no harm
  theoretical_peak_set.Clear();
  Peptide* peptide = getPeptide(i);
  peptide->ComputeTheoreticalPeaks(&theoretical_peak_set, dia_mode_);
}

bool PeptideDiskLoader::isWithinIsotope(vector<double>* min_mass, vector<double>* max_mass, double mass, int* isotope_idx) {
  for (int i = *isotope_idx; i < min_mass->size(); ++i) {
    if (mass >= (*min_mass)[i] && mass <= (*max_mass)[i]) {
      if (i > *isotope_idx) {
        *isotope_idx = i;
      }
      return true;
    }
  }
  return false;
}

const std::vector<RollingPeptideWindow*> PeptideDiskLoader::GetActivePeptideWindows() const {
  return windows_;
}

int RollingPeptideWindow::SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range) {
  //min_range and max_range have been introduced to fix a bug
  //introduced by m/z selection. see #222 in sourceforge
  //this has to be true:
  // min_range <= min_mass <= max_mass <= max_range

  // queue front() is lightest; back() is heaviest

  // delete anything already loaded that falls below min_range
  while (!empty() && front()->Mass() < min_range) {
    assert(PopFront());
  }

  /*
  while (!queue_.empty() && queue_.front()->Mass() < min_range) {
    Peptide* peptide = queue_.front();
    queue_.pop_front();
    delete peptide;
  }
  */
  nPeptides = 0;
  nCandPeptides = 0;
  CandPeptidesTarget = 0;
  CandPeptidesDecoy = 0;

  // Enqueue all peptides that are not yet queued but are lighter than
  // max_range. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.
  bool done = false;
  //Modified for tailor score calibration method by AKF

  while (size() < queue_->min_candidates_ && back()->Mass() <= max_range) {
    if (!(done = PushBack())) { // reader is done
      break;
    }
  }
  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!empty() || done);

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  if (empty()) {
    return 0;
  }

  nPeptides = 0;
  active_end_ = begin_;

  while (active_end_ < end_) {
    Peptide* peptide = queue_->getPeptide(active_end_);
    assert(peptide);

    if (peptide->Mass() >= min_mass->front()) {
      break;
    }
    peptide->active_ = false;
    active_end_++;
    nPeptides++;
  }
  active_begin_ = begin_;
  return nCandPeptides;
}

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
    return true;
  }
  return false;
}

bool PeptideDiskLoader::popFront(RollingPeptideWindow* window) {
  std::lock_guard<boost::shared_mutex> lock(m_);
  if (window->begin_ == window->end_) {
    return false;
  }
  Peptide* peptide = getPeptide(window->begin_);
  peptide->Drop();

  if (peptide->GetDropsCounter() == windows_.size()) {
    delete peptide;
    queue_.pop_front();
  }
  window->begin_++;
}

Peptide* PeptideDiskLoader::getPeptide(size_t i) {
  m_.lock_shared();
  Peptide* peptide = queue_.at(i - begin_);
  m_.unlock_shared();
  return peptide;
}

Peptide* PeptideDiskLoader::getComputedPeptide(size_t i) {
  ComputeTheoreticalPeak(i);
  return getPeptide(i);
}