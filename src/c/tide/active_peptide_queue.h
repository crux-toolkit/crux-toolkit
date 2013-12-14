// Benjamin Diament
//
// An ActivePeptideQueue is constructed with a file of peptides of
// non-decreasing neutral mass, and a correpsonding set of proteins. With
// successive call to SetActiveRange(min_mass, max_mass) the ActivePeptideQueue
// reads in peptides and initializes them. It also discards any peptides in 
// memory lighter than min_mass.
//
// The ActivePeptideQueue maintains in memory a window of peptides that fall
// within  a given mass range. The window to maintain (the "active peptides")
// is set by a call to SetActiveRange(). Successive values
// of min_mass and max_mass must be non-decreasing. After the call to
// SetActiveRange() the client may use the iterator interface HasNext() and
// NextPeptide() to iterate over the window. The client may also use
// GetPeptide() to get a specific peptide in the window.

#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"

class TheoreticalPeakCompiler;

class ActivePeptideQueue {
 public:
  ActivePeptideQueue(RecordReader* reader,
		     const vector<const pb::Protein*>& proteins);

  ~ActivePeptideQueue();

  // See above for usage and .cc for implmentation details.
  int SetActiveRange(double min_mass, double max_mass);

  bool HasNext() const { return iter_ != end_; }
  Peptide* NextPeptide() { return *iter_++; }
  const Peptide* GetPeptide(int back_index) const { return *(end_ - back_index); }

  int ActiveTargets() const { return active_targets_; }
  int ActiveDecoys() const { return active_decoys_; }

 // IMPLEMENTATION DETAILS
 
 private:
  // See .cc file.
  void ComputeTheoreticalPeaksBack();

  RecordReader* reader_;
  pb::Peptide current_pb_peptide_;
  
  // All amino acid sequences from which the peptides are drawn.
  const vector<const pb::Protein*>& proteins_; 

  // Workspace for computing theoretical peaks for a single peptide.
  // Gets reused for each new peptide.
  ST_TheoreticalPeakSet theoretical_peak_set_;
  
  // The active peptides. Lighter peptides are enqueued before heavy ones.
  // queue_ maintains only the peptides that fall within the range specified
  // by the last call to SetActiveRange().
  deque<Peptide*> queue_;

  // iter_ points to the current peptide. Client access is by HasNext(),
  // GetPeptide(), and NextPeptide(). end_ points just beyond the last active
  // peptide.
  deque<Peptide*>::const_iterator iter_, end_;

  // Set by most recent call to SetActiveRange()
  double min_mass_, max_mass_;

  // While we maintain a window of active peptides, we allocate and relase them
  // on a first-in, first-out basis. We use FifoAllocators 
  // (see fifo_alloc.{h,cc}) to manage memory efficiently for this usage 
  // pattern. As we read in peptides and compute the theoretical peaks, we
  // use compiler_prog1 and compiler_prog2 to generate code on the fly for
  // taking dot products with the theoretical peak set for each peptide. 
  // FifoAllocators allow us to execute the code thus generated,
  // since they set the proper permissions. The set of theoretical peaks for 
  // "dotting" with charge 1 and charge 2 spectra, have different
  // FifoAllocators and TheoreticalPeakCompilers.
  FifoAllocator fifo_alloc_peptides_;
  FifoAllocator fifo_alloc_prog1_;
  FifoAllocator fifo_alloc_prog2_;
  TheoreticalPeakCompiler* compiler_prog1_;
  TheoreticalPeakCompiler* compiler_prog2_;

  // Number of targets and decoys in active range
  int active_targets_, active_decoys_;
};

/*
CONSIDER:
Make 3 subclasses:
  Read-as-you-go (below)
  Preread into memory (as current operation)
  Threaded reads (!!)

pb::Peptide* pb_peptide = new pb::Peptide;
CHECK(reader_->Read(pb_peptide));
!reader_->Done()
*/
