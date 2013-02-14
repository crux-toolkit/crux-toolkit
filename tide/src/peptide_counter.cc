#include <iostream>
#include <gflags/gflags.h>
#include "records.h"
#include "records_to_vector-inl.h"
#include "peptides.pb.h"

using namespace std;

DEFINE_string(peptides, "", "File of peptides to count, as peptides.proto");

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  int count = CountRecords<pb::Peptide>(FLAGS_peptides);

  if (count < 0)
    return count;

  cout << count << " RECORDS" << endl;

  return 0;
}
