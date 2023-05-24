// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
#include <deque>
#include <gflags/gflags.h>
#include "records.h"
#include "peptides.pb.h"
#include "peptide.h"
#include "ActivePeptideQueueLite.h"
#include "records_to_vector-inl.h"
#include "theoretical_peak_set.h"
#include "compiler.h"
#include "app/TideLiteMatchSet.h"
#include <map> //Added by Andy Lin
#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_int32(fifo_page_size, 1, "Page size for FIFO allocator, in megs");

ActivePeptideQueueLite::ActivePeptideQueueLite(RecordReader* reader,
                                       const vector<const pb::Protein*>& proteins, 
                                       vector<const pb::AuxLocation*>* locations) {
}

ActivePeptideQueueLite::~ActivePeptideQueueLite() {
}
