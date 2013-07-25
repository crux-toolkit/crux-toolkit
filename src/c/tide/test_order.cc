#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <gflags/gflags.h>
#include "records.h"
#include "header.pb.h"
#include "peptides.pb.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

#if 0

DEFINE_string(in, "", "File to read");

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  pb::Peptide p;

  RecordReader r(FLAGS_in);
  CHECK(r.OK());

  double last = 0;
  for (int i = 0; !r.Done(); ++i) {
    r.Read(&p);
    double m = p.mass();
    CHECK(last <= m) << "INVERSION AT " << i;
    last = m;
    if (i % 1000000 == 0)
      fprintf(stderr, "%d\r", i);
  }
  printf("\n");

  CHECK(r.OK());

  return 0;
}

#endif

DEFINE_string(a, "", "First file");
DEFINE_string(b, "", "Second file");

bool IsNear(double x, double y) {
  const double tol = 1e-9;
  x-=y;
  return x<tol && x>-tol;
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  pb::Peptide pa, pb;

  RecordReader ra(FLAGS_a, (1<<20));
  RecordReader rb(FLAGS_b, (1<<20) + 170000);
  CHECK(ra.OK());
  CHECK(rb.OK());

  for (int i = 0; !ra.Done(); ++i) {
    CHECK(!rb.Done());
    ra.Read(&pa);
    rb.Read(&pb);
    CHECK(pa.id() == pb.id());
    CHECK(pa.id() == pb.id());
    CHECK(pa.length() == pb.length());
    CHECK(pa.first_location().protein_id() == pb.first_location().protein_id());
    CHECK(pa.first_location().pos() == pb.first_location().pos());
    CHECK(pa.modifications_size() == pb.modifications_size());
    for (int j = 0; j < pa.modifications_size(); ++j)
      CHECK(pa.modifications(j) == pb.modifications(j));
    //cout << pa.mass() << "   " <<  pb.mass() << "   " << (pa.mass() - pb.mass()) << endl;
    CHECK(IsNear(pa.mass(), pb.mass()));
    if (i % 1000000 == 0)
      fprintf(stderr, "%d\r", i);
  }
  printf("\n");

  CHECK(ra.OK());
  CHECK(rb.OK());

  return 0;
}

