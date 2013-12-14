// Benjamin Diament
#include <deque>
#include <gflags/gflags.h>
#include "records.h"
#include "peptides.pb.h"
#include "peptide.h"
#include "active_peptide_queue.h"
#include "records_to_vector-inl.h"
#include "theoretical_peak_set.h"
#include "compiler.h"

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_int32(fifo_page_size, 1, "Page size for FIFO allocator, in megs");

ActivePeptideQueue::ActivePeptideQueue(RecordReader* reader,
                                       const vector<const pb::Protein*>&
                                       proteins)
  : reader_(reader),
    proteins_(proteins),
    theoretical_peak_set_(2000), // probably overkill, but no harm
    active_targets_(0), active_decoys_(0),
    fifo_alloc_peptides_(FLAGS_fifo_page_size << 20),
    fifo_alloc_prog1_(FLAGS_fifo_page_size << 20),
    fifo_alloc_prog2_(FLAGS_fifo_page_size << 20) {
  CHECK(reader_->OK());
  compiler_prog1_ = new TheoreticalPeakCompiler(&fifo_alloc_prog1_);
  compiler_prog2_ = new TheoreticalPeakCompiler(&fifo_alloc_prog2_);
}

ActivePeptideQueue::~ActivePeptideQueue() {
  deque<Peptide*>::iterator i = queue_.begin();
  // for (; i != queue_.end(); ++i)
  //   delete (*i)->PB();
  fifo_alloc_peptides_.ReleaseAll();
  fifo_alloc_prog1_.ReleaseAll();
  fifo_alloc_prog2_.ReleaseAll();

  delete compiler_prog1_;
  delete compiler_prog2_;
}

// Compute the theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueue::ComputeTheoreticalPeaksBack() {
  theoretical_peak_set_.Clear();
  Peptide* peptide = queue_.back();
  peptide->ComputeTheoreticalPeaks(&theoretical_peak_set_, current_pb_peptide_,
                                   compiler_prog1_, compiler_prog2_);
}

int ActivePeptideQueue::SetActiveRange(double min_mass, double max_mass) {
  // queue front() is lightest; back() is heaviest
  
  // delete anything already loaded that falls below min_mass
  while (!queue_.empty() && queue_.front()->Mass() < min_mass) {
    Peptide* peptide = queue_.front();
    // would delete peptide's underlying pb::Peptide;
    queue_.pop_front();
  }
  if (queue_.empty()) {
    //cerr << "Releasing All\n";
    fifo_alloc_peptides_.ReleaseAll();
    fifo_alloc_prog1_.ReleaseAll();
    fifo_alloc_prog2_.ReleaseAll();
    //cerr << "Prog1: ";
    //fifo_alloc_prog1_.Show();
    //cerr << "Prog2: ";
    //fifo_alloc_prog2_.Show();
  } else {
    Peptide* peptide = queue_.front();
    // Free all peptides up to, but not including peptide.
    fifo_alloc_peptides_.Release(peptide); 
    peptide->ReleaseFifo(&fifo_alloc_prog1_, &fifo_alloc_prog2_);
  }
  
  // Enqueue all peptides that are not yet queued but are lighter than
  // max_mass. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.
  bool done = false;
  if (queue_.empty() || queue_.back()->Mass() <= max_mass) {
    if (!queue_.empty())
      ComputeTheoreticalPeaksBack();
    while (!(done = reader_->Done())) {
      // read all peptides lighter than max_mass
      reader_->Read(&current_pb_peptide_);
      if (current_pb_peptide_.mass() < min_mass) {
        // we would delete current_pb_peptide_;
        continue; // skip peptides that fall below min_mass
      }
      Peptide* peptide = new(&fifo_alloc_peptides_)
        Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);
      queue_.push_back(peptide);
      if (peptide->Mass() > max_mass)
        break;
      ComputeTheoreticalPeaksBack();
    }
  }
  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!queue_.empty() || done);
  
  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  iter_ = queue_.begin();
  end_ = queue_.end();
  if (queue_.empty())
    return 0;

  int active = queue_.size();
  if (!done) {
    --end_;
    --active;
  }
  // Count active targets and decoys
  active_targets_ = active_decoys_ = 0;
  for (deque<Peptide*>::const_iterator i = iter_; i != end_; ++i) {
    if (!(*i)->IsDecoy()) {
      ++active_targets_;
    } else {
      ++active_decoys_;
    }
  }
  return active;

  /*
  cerr << (end_ - iter_) << " candidates.";
  if (end_ != iter_)
    cerr << " Range: (" << (*iter_)->PB()->id() << ", " << (*end_)->PB()->id() << ")";
  cerr << endl;
  */
}
