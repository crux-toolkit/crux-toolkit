// Benjamin Diament
//
// An indexing phase to sort a file of Peptides by mass. To be done after
// make_peptides.py (tryptic digestion). 
//
// Example command-line:
// peptide_sorter --in_peptides=<peptides.proto input file> \
//                --out_peptides=<peptides.proto output file>
//
// Peptides are all read into memory and then sorted; not scalable (TODO 249).

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "records.h"
#include "records_to_vector-inl.h"
#include "peptides.pb.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(in_peptides, "", "File of unsorted peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to create of sorted peptides, as peptides.proto");

struct less_peptide : public binary_function<const pb::Peptide*, const pb::Peptide*, bool> {
  bool operator()(const pb::Peptide* x, const pb::Peptide* y) { return x->mass() < y->mass(); }
};

void SortPeptides(const string& input_filename,
		  const string& output_filename) {
  RecordWriter writer(output_filename);
  CHECK(writer.OK());

  typedef vector<pb::Peptide*> Vec;
  Vec peptides;

  CHECK(ReadRecordsToVector<pb::Peptide>(&peptides, input_filename));
  
  sort(peptides.begin(), peptides.end(), less_peptide());

  for (int i = 0; i < peptides.size(); ++i) {
    peptides[i]->set_id(i);
    CHECK(writer.Write(peptides[i]));
  }
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  SortPeptides(FLAGS_in_peptides, FLAGS_out_peptides);

  return 0;
}
