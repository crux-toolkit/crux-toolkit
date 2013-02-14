// Benjamin Diament
//
// Add to the index of peptide records the pre-computed theoretical peaks.
// We store the TheoreticalPeakSetDiff (q.v.) for each peptide. 
//
// Example command-line:
// peptide_peaks --proteins=<raw_proteins.proto input file> \
//               --in_peptides=<peptides.proto input file> \
//               --out_peptides=<peptides.proto output file> \

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <gflags/gflags.h>
#include "peptides.pb.h"
#include "records.h"
#include "records_to_vector-inl.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

DEFINE_string(proteins, "", "File of proteins corresponding to peptides, as raw_proteins.proto");
DEFINE_string(in_peptides, "", "File of unfragmented peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to create of fragmented peptides, as peptides.proto");

extern void AddTheoreticalPeaks(const vector<const pb::Protein*>& proteins,
				const string& input_filename,
				const string& output_filename);

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  vector<const pb::Protein*> proteins;
  CHECK(ReadRecordsToVector<pb::Protein>(&proteins, FLAGS_proteins));
  AddTheoreticalPeaks(proteins, FLAGS_in_peptides, FLAGS_out_peptides);

  return 0;
}
