#include <iostream>
#include <fstream>
#include <vector>
#include <gflags/gflags.h>
#include "spectrum.pb.h"
#include "spectrum_collection.h"
#include "records.h"
#include "abspath.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK((x))

DEFINE_string(in, "", "Input spectra");
DEFINE_string(out, "", "Filtered output spectra");
DEFINE_int32(min_peaks, 5, "Filter out spectra with fewer peaks");

bool Filter(const Spectrum& spec) {
  return spec.Size() >= FLAGS_min_peaks;
}

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  pb::Header header;
  HeadedRecordReader reader(FLAGS_in, &header);
  CHECK(header.file_type() == pb::Header::SPECTRA);

  pb::Header new_header;
  new_header.set_file_type(pb::Header::SPECTRA);
  pb::Header_Source* source = new_header.add_source();
  source->set_filename(AbsPath(FLAGS_in));
  source->mutable_header()->CopyFrom(header);
  HeadedRecordWriter writer(FLAGS_out, new_header);
  

  int in_count = 0;
  int out_count = 0;
  pb::Spectrum pb_spectrum;
  while (!reader.Done()) {
    reader.Read(&pb_spectrum);
    ++in_count;
    Spectrum spec(pb_spectrum);
    if (Filter(spec)) {
      CHECK(writer.Write(&pb_spectrum));
      ++out_count;
    }
  }
  CHECK(reader.OK());
  CHECK(writer.OK());

  cout << in_count << " input spectra.\n"
       << out_count << " filtered output spectra.\n";

  return 0;
}
