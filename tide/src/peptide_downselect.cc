// Benjamin Diament
//
// Randomly sample from a collection of peptides. See command line flags
// below.
// Both input and output files are records of protocol buffer messages as
// peptides.proto (q.v.).

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "records.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(in_peptides, "", "File of sorted peptides, as peptides.proto");
DEFINE_string(out_peptides, "", "File to create of randomly-selected peptides, as peptides.proto");
DEFINE_double(fraction, 0.02, "Fraction of peptides to select");

inline bool Select() {
  static int threshold = int(RAND_MAX * FLAGS_fraction + 0.5);
  return random() < threshold;
}

void SelectPeptides(const string& input_filename,
		    const string& output_filename) {
  pb::Header header;
  HeadedRecordReader reader(input_filename, &header);
  CHECK(header.file_type() == pb::Header::PEPTIDES);
  CHECK(header.has_peptides_header());
  header.mutable_peptides_header()->set_downselect_fraction(FLAGS_fraction);

  HeadedRecordWriter writer(output_filename, header);
  CHECK(writer.OK());

  pb::Peptide peptide;
  int count = 0;
  while (!reader.Done()) {
    reader.Read(&peptide);
    if (++count % 100000 == 0) // Report progress once in a while.
      cerr << "Read " << count << " peptides\n";
    if (Select())
      CHECK(writer.Write(&peptide));
  }
  CHECK(reader.OK());
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  SelectPeptides(FLAGS_in_peptides, FLAGS_out_peptides);

  return 0;
}
