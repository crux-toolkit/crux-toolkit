#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <gflags/gflags.h>
#include "spectrum.pb.h"
#include "spectrum_collection.h"
#include "records.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(out, "", "Combined output spectra");

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  int inputs[] = {90,93,94,95,96,97};
  int num_inps = sizeof(inputs)/sizeof(inputs[0]);
  pb::Header header[num_inps];
  HeadedRecordReader* reader[num_inps];
  for (int inp = 0; inp < num_inps; ++inp) {
    stringstream name;
    name << inputs[inp] << ".spectrumrecords";
    reader[inp] = new HeadedRecordReader(name.str(), &header[inp]);
    CHECK(header[inp].file_type() == pb::Header::SPECTRA);
  }

  pb::Header new_header;
  new_header.set_file_type(pb::Header::SPECTRA);
  for (int inp = 0; inp < num_inps; ++inp) {
    pb::Header_Source* source = new_header.add_source();
    stringstream name;
    name << inputs[inp] << ".spectrumrecords";
    source->set_filename(AbsPath(name.str()));
    source->mutable_header()->CopyFrom(header[inp]);
  }
  HeadedRecordWriter writer(FLAGS_out, new_header);

  int count = 0;
  pb::Spectrum pb_spectrum;
  for (int inp = 0; inp < num_inps; ++inp) {
    while (!reader[inp]->Done()) {
      reader[inp]->Read(&pb_spectrum);
      ++count;
      int specnum = inputs[inp] * 10000000 + pb_spectrum.spectrum_number();
      pb_spectrum.set_spectrum_number(specnum);
      CHECK(writer.Write(&pb_spectrum));
    }
    CHECK(reader[inp]->OK());
  }
  CHECK(writer.OK());

  cout << count << " spectra.\n";

  return 0;
}
