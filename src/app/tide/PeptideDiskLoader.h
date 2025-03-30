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

#ifndef ACTIVE_PEPTIDE_QUEUE_H
#define ACTIVE_PEPTIDE_QUEUE_H

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

  const std::vector<RollingPeptideWindow*> GetActivePeptideWindows() const;

  std::deque<Peptide*> queue_;
  int min_candidates_;
  bool dia_mode_;

 private:
  bool popFront(RollingPeptideWindow*);
  bool pushBack(RollingPeptideWindow*);


  static bool isWithinIsotope(vector<double>* min_mass,
        vector<double>* max_mass, 
        double mass,
        int* isotope_idx);   

  void ComputeTheoreticalPeak(size_t i);

  Peptide* getPeptide(size_t i);

  Peptide* getComputedPeptide(size_t i);
  
  boost::shared_mutex m_;

  size_t begin_ = 1;
  size_t end_ = 1;


  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_; 
  vector<const pb::AuxLocation*>* locations_;
  
  pb::Peptide current_pb_peptide_;

  std::vector<RollingPeptideWindow*> windows_;
};

class RollingPeptideWindow {
public:
  friend class PeptideDiskLoader;

  int SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, 
      double min_range, double max_range);
  bool PushBack();
  bool PopFront();

  inline Peptide* back() { return (begin_ == end_) ? NULL : queue_->getPeptide(end_ - 1); }
  inline Peptide* front() { return (begin_ == end_) ? NULL : queue_->getPeptide(begin_); }

  inline bool empty() const { return begin_ == end_; }
  inline size_t size() const { return end_ - begin_; }

  inline size_t begin() { return begin_; }
  inline size_t end() { return end_; }

  inline size_t active_begin() { return active_begin_; }
  inline size_t active_end() { return active_end_;}

  Peptide* GetPeptide(size_t i) { return queue_->getPeptide(i); }

  PeptideDiskLoader* GetQueue() const { return queue_; }


  int nPeptides= 0;
  int nCandPeptides = 0;
  int CandPeptidesTarget = 0;
  int CandPeptidesDecoy = 0;
private:
  RollingPeptideWindow(PeptideDiskLoader* queue);

  size_t begin_ = 1;
  size_t end_ = 1;
  size_t active_begin_ = 0;
  size_t active_end_ = 0;
  PeptideDiskLoader* queue_;
};
#endif
