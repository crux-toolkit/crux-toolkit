// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
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
    theoretical_peak_set_(2000),   // probably overkill, but no harm
    theoretical_b_peak_set_(200),  // probably overkill, but no harm
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

// Compute the b ion only theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueue::ComputeBTheoreticalPeaksBack() {
  theoretical_b_peak_set_.Clear();
  Peptide* peptide = queue_.back();
  peptide->ComputeBTheoreticalPeaks(&theoretical_b_peak_set_);
  b_ion_queue_.push_back(theoretical_b_peak_set_);
}

int ActivePeptideQueue::SetActiveRangeBIons(double min_mass, double max_mass) {
  // queue front() is lightest; back() is heaviest
  
  // delete anything already loaded that falls below min_mass
  while (!queue_.empty() && queue_.front()->Mass() < min_mass) {
    Peptide* peptide = queue_.front();
    // would delete peptide's underlying pb::Peptide;
    queue_.pop_front();
    b_ion_queue_.pop_front();
  }
  if (queue_.empty()) {
    fifo_alloc_peptides_.ReleaseAll();
  } else {
    Peptide* peptide = queue_.front();
    // Free all peptides up to, but not including peptide.
    fifo_alloc_peptides_.Release(peptide); 
  }
  
  // Enqueue all peptides that are not yet queued but are lighter than
  // max_mass. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.
  bool done;
  if (queue_.empty() || queue_.back()->Mass() <= max_mass) {
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
      ComputeBTheoreticalPeaksBack();
      if (peptide->Mass() > max_mass)
        break;
    }
  }
  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!queue_.empty() || done);
  
  // Set up iterators for use with b_ion_queue_
  iter1_ = b_ion_queue_.begin();
  end1_ = b_ion_queue_.end();
  --end1_;

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  iter_ = queue_.begin();
  end_ = queue_.end();
  if (queue_.empty())
    return 0;
  --end_;

  // Count active targets and decoys
  active_targets_ = active_decoys_ = 0;
  for (deque<Peptide*>::const_iterator i = iter_; i != end_; ++i) {
    if (!(*i)->IsDecoy()) {
      ++active_targets_;
    } else {
      ++active_decoys_;
    }
  }

  int active = queue_.size();
  --active;
  return active;

  /*
  cerr << (end_ - iter_) << " candidates.";
  if (end_ != iter_)
    cerr << " Range: (" << (*iter_)->PB()->id() << ", " << (*end_)->PB()->id() << ")";
  cerr << endl;
  */
}

int ActivePeptideQueue::CountAAFrequency(
  double binWidth,
  double binOffset,
  double** dAAFreqN,
  double** dAAFreqI,
  double** dAAFreqC,
  int** dAAMass
) {

    unsigned int i = 0;
    unsigned int cntTerm = 0;
    unsigned int cntInside = 0;
    const unsigned int MaxModifiedAAMassBin = 2000 / binWidth;   //2000 is the maximum size of a modified amino acid 
    unsigned int* nvAAMassCounterN = new unsigned int[MaxModifiedAAMassBin];   //N-terminal amino acids
    unsigned int* nvAAMassCounterC = new unsigned int[MaxModifiedAAMassBin];   //C-terminal amino acids
    unsigned int* nvAAMassCounterI = new unsigned int[MaxModifiedAAMassBin];   //inner amino acids in the peptides
    memset(nvAAMassCounterN, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
    memset(nvAAMassCounterC, 0, MaxModifiedAAMassBin * sizeof(unsigned int));
    memset(nvAAMassCounterI, 0, MaxModifiedAAMassBin * sizeof(unsigned int));

    while (!(reader_->Done())) { // read all peptides in index
      reader_->Read(&current_pb_peptide_);
      Peptide* peptide = new(&fifo_alloc_peptides_) Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);

      double* dAAResidueMass = peptide->getAAMasses(); //retrieves the amino acid masses, modifications included

      int nLen = peptide->Len(); //peptide length
      ++nvAAMassCounterN[(unsigned int)(dAAResidueMass[0] / binWidth + 1.0 - binOffset)];
      for (i = 1; i < nLen-1; ++i) {  
        ++nvAAMassCounterI[(unsigned int)(dAAResidueMass[i] / binWidth + 1.0 - binOffset)];
        ++cntInside;
      } 
      ++nvAAMassCounterC[(unsigned int)(dAAResidueMass[nLen - 1] / binWidth + 1.0 - binOffset)];
      ++cntTerm;

      delete dAAResidueMass;
      fifo_alloc_peptides_.ReleaseAll();
    }

  //calculate the unique masses
  unsigned int uiUniqueMasses = 0;
  for (i = 0; i < MaxModifiedAAMassBin; ++i) {
    if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
      ++uiUniqueMasses;
    }
  }

  //calculate the unique amino acid masses
  *dAAMass = new int[uiUniqueMasses];     //a vector for the unique (integerized) amino acid masses present in the sample
  *dAAFreqN = new double[uiUniqueMasses]; //a vector for the amino acid frequencies at the N-terminus
  *dAAFreqI = new double[uiUniqueMasses]; //a vector for the amino acid frequencies inside the peptide
  *dAAFreqC = new double[uiUniqueMasses]; //a vector for the amino acid frequencies at the C-terminus
  unsigned int cnt = 0;
  for (i = 0; i < MaxModifiedAAMassBin; ++i) {
    if (nvAAMassCounterN[i] || nvAAMassCounterI[i] || nvAAMassCounterC[i]) {
      (*dAAFreqN)[cnt] = (double)nvAAMassCounterN[i] / cntTerm;
      (*dAAFreqI)[cnt] = (double)nvAAMassCounterI[i] / cntInside;
      (*dAAFreqC)[cnt] = (double)nvAAMassCounterC[i] / cntTerm;
      (*dAAMass)[cnt] = i;
      cnt++;
    }
  }

  delete nvAAMassCounterN;
  delete nvAAMassCounterI;
  delete nvAAMassCounterC;
  return uiUniqueMasses;
}
