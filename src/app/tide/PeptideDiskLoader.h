#include <deque>
#include <list>
#include <unordered_set>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"
#include "boost/thread/shared_mutex.hpp"
#include "boost/thread/shared_lock_guard.hpp"

#ifndef PEPTIDE_DISK_LOADER_H
#define PEPTIDE_DISK_LOADER_H

class TheoreticalPeakCompiler;
class RollingPeptideWindow;

class PeptideDiskLoader {
 public:
  using iterator = std::deque<Peptide*>::iterator;
  using const_iterator = std::deque<Peptide*>::const_iterator;
  friend class RollingPeptideWindow;

  PeptideDiskLoader(RecordReader* reader,
        const vector<const pb::Protein*>& proteins,
        vector<const pb::AuxLocation*>* locations=NULL, 
        bool dia_mode = false, int thread_num = 1);

  ~PeptideDiskLoader();

  const std::vector<RollingPeptideWindow*> GetRollingPeptideWindows() const;

  std::deque<Peptide*> queue_;
  int min_candidates_;
  bool dia_mode_;

 private:
  bool popFront(RollingPeptideWindow*);
  bool pushBack(RollingPeptideWindow*);

  // void ComputeTheoreticalPeak(size_t i);
  // void ComputeTheoreticalPeak(Peptide* peptide);

  Peptide* getPeptide(size_t i);
  Peptide* getPeptideUnsafe(size_t i);

  // Peptide* getComputedPeptide(size_t i);
  
  boost::shared_mutex m_;

  size_t end_ = 0;
  size_t begin_ = 0;

  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_; 
  vector<const pb::AuxLocation*>* locations_;
  
  pb::Peptide current_pb_peptide_;

  std::vector<RollingPeptideWindow*> windows_;
};

class RollingPeptideWindow {
public:
  friend class PeptideDiskLoader;

  int SetActiveRange(double min_range, double max_range, vector<double>* min_mass, vector<double>* max_mass);
  bool PushBack();
  bool PopFront();

  inline Peptide* back() { return (begin_ == end_) ? NULL : queue_->getPeptide(end_ - 1); }
  inline Peptide* front() { return (begin_ == end_) ? NULL : queue_->getPeptide(begin_); }

  inline bool empty() const { return size() == 0; }
  inline size_t size() const { return end_ - begin_; }

  inline size_t begin() { return begin_; }
  inline size_t end() { return end_; }

  Peptide* GetPeptide(size_t i) { return queue_->getPeptide(i); }

  PeptideDiskLoader* GetQueue() const { return queue_; }

  size_t nPeptides_ = 0;
  size_t nCandPeptides_ = 0;
  size_t CandPeptidesTarget_ = 0;
  size_t CandPeptidesDecoy_ = 0;

  vector<double>* min_mass_;
  vector<double>* max_mass_;

  TheoreticalPeakSetBYSparse theoretical_peak_set_; // Probably overkill, but no harm    
private:
  RollingPeptideWindow(PeptideDiskLoader* queue);

  size_t begin_ = 0;
  size_t end_ = 0;
  PeptideDiskLoader* queue_;
};
#endif
