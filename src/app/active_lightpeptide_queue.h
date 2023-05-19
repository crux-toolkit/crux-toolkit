#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"

#ifndef ACTIVE_PEPTIDE_QUEUE_H
#define ACTIVE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

class ActivePeptideQueue {
 public:
  ActivePeptideQueue(); //define in active_peptide_queue.h
  ~ActivePeptideQueue();
  
 private:
};
#endif //ACTIVE_PEPTIDE_QUEUE_H
