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

DEFINE_string(in, "", "MS2 File of spectra to read");
DEFINE_string(out, "", "Protocol Buffer file of spectra to write");
DEFINE_bool(ms1, false, "Parse as ms1 instead of ms2");
DEFINE_bool(sort, false, "Sort spectra before writing. Will write each "
	    "spectrum once for each charge state, so output will be larger "
	    "than input");

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  ifstream in(FLAGS_in.c_str());
  SpectrumCollection spectrum_collection;
  spectrum_collection.ReadMS(in, FLAGS_ms1);
  in.close();
  vector<Spectrum*>* spectra = spectrum_collection.Spectra();

  cout << spectra->size() << " SPECTRA" << endl;

  if (!FLAGS_out.empty()) {
    pb::Header header;
    header.set_file_type(pb::Header::SPECTRA);
    pb::Header_Source* source = header.add_source();
    source->set_filename(AbsPath(FLAGS_in));
    source->set_filetype("ms2");
    header.mutable_spectra_header()->set_sorted(FLAGS_sort);
    HeadedRecordWriter writer(FLAGS_out, header);
    CHECK(writer.OK());
    if (FLAGS_sort) {
      spectrum_collection.Sort();
      const vector<SpectrumCollection::SpecCharge>* spec_charges
	= spectrum_collection.SpecCharges();
      for (int i = 0; i < spec_charges->size(); ++i) {
	pb::Spectrum spectrum;
	(*spec_charges)[i].spectrum->FillPB(&spectrum);
	spectrum.clear_charge_state();
	spectrum.add_charge_state((*spec_charges)[i].charge);
	CHECK(writer.Write(&spectrum));
      }
    } else {
      pb::Spectrum spectrum;
      for (int i = 0; i < spectra->size(); ++i) {
	(*spectra)[i]->FillPB(&spectrum);
	CHECK(writer.Write(&spectrum));
      }
    }
  }

  return 0;
}
