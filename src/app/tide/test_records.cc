// Benjamin Diament
//
// Write and/or read a set of 100 made-up raw_proteins.proto records.
// The contents exactly match those tested in test-records.py.

#include <stdio.h>
#include <vector>
#include <gflags/gflags.h>
#include "header.pb.h"
#include "raw_proteins.pb.h"
#include "records.h"
#include "records_to_vector-inl.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(testfile, "/tmp/testrecords", "temp file to create and read");

pb::Header* GetTestHeader() {
  pb::Header* header = new pb::Header();
  header->set_file_type(pb::Header::RAW_PROTEINS);
  pb::Header_Source* source = header->add_source();
  source->set_filename("this would be the name of a fasta file");
  source->set_filetype("fasta");
  return header;
}

vector<pb::Protein*>* GetTestProteins() {
  vector<pb::Protein*>* vec = new vector<pb::Protein*>;
  for (int id = 0; id < 100; ++id) {
    pb::Protein* protein = new pb::Protein;
    protein->set_id(id);
    protein->set_name("some name");
    protein->set_residues("RESID" "LETTERAFTERT" "ES");
    vec->push_back(protein);
  }
  return vec;
}

void DeleteTestProteins(vector<pb::Protein*>* vec) {
  for (vector<pb::Protein*>::iterator i = vec->begin(); i != vec->end(); ++i)
    delete (*i);
}

void TestWrite(const pb::Header* header, const vector<pb::Protein*>* vec) {
  HeadedRecordWriter writer(FLAGS_testfile, *header);
  CHECK(writer.OK());

  vector<pb::Protein*>::const_iterator i = vec->begin();
  for (; i != vec->end(); ++i)
    writer.Write(*i);

  CHECK(writer.OK());
}

#define EQ(x, y) { if ((x) != (y)) return false; }

bool Match(const pb::Protein& a, const pb::Protein& b) {
  EQ(a.id(), b.id());
  EQ(a.name(), b.name());
  EQ(a.residues(), b.residues());
  return true;
}

bool Match(const pb::Header& a, const pb::Header& b) {
  EQ(a.file_type(), b.file_type());
  EQ(a.source_size(), 1);
  EQ(b.source_size(), 1);
  EQ(a.source(0).filename(), b.source(0).filename());
  EQ(a.source(0).filetype(), b.source(0).filetype());
  return true;
}

void TestRead(const pb::Header* correct_header,
	      const vector<pb::Protein*>* correct_vec) {
  pb::Header test_header;
  vector<pb::Protein*> test_vec;
  CHECK(ReadRecordsToVector<pb::Protein>(&test_vec, FLAGS_testfile,
					 &test_header));
  CHECK(Match(test_header, *correct_header));
  CHECK(test_vec.size() == correct_vec->size());
  vector<pb::Protein*>::const_iterator i = test_vec.begin();
  vector<pb::Protein*>::const_iterator j = correct_vec->begin();
  for (; i != test_vec.end(); ++i, ++j)
    CHECK(Match(**i, **j));
  DeleteTestProteins(&test_vec);
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  pb::Header* header = GetTestHeader();
  vector<pb::Protein*>* vec = GetTestProteins();
  
  TestWrite(header, vec);
  TestRead(header, vec);

  DeleteTestProteins(vec);
  delete header;
  delete vec;

  printf("PASSED.\n");

  return 0;
}
