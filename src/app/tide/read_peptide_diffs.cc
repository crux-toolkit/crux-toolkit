#include <stdio.h>
#include <gflags/gflags.h>
#include "records.h"
#include "raw_proteins.pb.h"
#include "peptides.pb.h"
#include "theoretical_peak_pair.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(peptides, "", "File of unfragmented peptides, as "
	      "peptides.proto");


int g_threes = 0;
int g_all = 0;
int g_diffs[100000] = {0};

void Show(const google::protobuf::RepeatedField<int>& src) {
  int total = 0;
  google::protobuf::RepeatedField<int>::const_iterator i = src.begin();
  for (; i != src.end(); ++i) {
    if (total != 0) {
      assert(*i < 100000);
      ++g_diffs[*i];
    }
    TheoreticalPeakPair peak(total += *i);
    printf(" (%d, %d)", peak.Bin(), peak.Type());
    if (peak.Type() == 3)
      ++g_threes;
  }
  g_all += src.size();
}

void ReadPeptides(const string& filename) {
  RecordReader reader(filename);
  pb::Peptide peptide;
  while (!reader.Done()) {
    reader.Read(&peptide);
    printf("Peptide ID: %d  len: %d\n", peptide.id(), peptide.length());

    printf("  peak1:");
    Show(peptide.peak1());
    printf("\n");

    printf("  peak2:");
    Show(peptide.peak2());
    printf("\n");

    printf("  neg_peak1:");
    Show(peptide.neg_peak1());
    printf("\n");

    printf("  neg_peak2:");
    Show(peptide.neg_peak2());
    printf("\n");

    printf("\n");
  }
  CHECK(reader.OK());
}


int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);
  
  CHECK(!FLAGS_peptides.empty());

  ReadPeptides(FLAGS_peptides);

  printf("%d/%d threes\n", g_threes, g_all);
  int num_diffs = sizeof(g_diffs) / sizeof(g_diffs[0]);
  for (int i = 0; i < num_diffs; ++i)
    if (g_diffs[i] > 0)
      printf("Diff %d Count %d\n", i, g_diffs[i]);

  return 0;
}
