#include <deque>
#include "peptides.pb.h"
#include "peptide.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"

#ifndef ACTIVE_LITE_PEPTIDE_QUEUE_H
#define ACTIVE_LITE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

class ActivePeptideQueueLite {
 public:
    ActivePeptideQueueLite(RecordReader* reader,
            const vector<const pb::Protein*>& proteins,
            vector<const pb::AuxLocation*>* locations=NULL);
     ~ActivePeptideQueueLite();
 private:
};

#endif
