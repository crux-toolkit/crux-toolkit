#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <gflags/gflags.h>
#include "spectrum_collection.h"
#include "records.h"

using namespace std;

void Show(HeadedRecordReader& reader) {
  int count = 0;
  pb::Spectrum pb_spectrum;
  cout << setprecision(10);
  while (!reader.Done()) {
    reader.Read(&pb_spectrum);
    Spectrum spectrum(pb_spectrum);
    
    cout << "Spectrum Number: " << spectrum.SpectrumNumber();
    if (spectrum.PrecursorMZ() > 0)
      cout << "  Precursor m/z: " << spectrum.PrecursorMZ();
    cout << endl;

    if (spectrum.NumChargeStates() > 0) {
      cout << "Charge states: " << spectrum.ChargeState(0);
      for (int i = 1; i < spectrum.NumChargeStates(); ++i)
	cout << ", " << spectrum.ChargeState(i);
      cout << endl;
    }

    for (int i = 0; i < spectrum.Size(); ++i)
      cout << spectrum.M_Z(i) << " " << spectrum.Intensity(i) << endl;
  }
}

#define CHECK(x) GOOGLE_CHECK(x)

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (argc != 2) {
    cerr << "Usage:  " << argv[0] << " SPECTRUM_FILE" << endl;
    return -1;
  }

  pb::Header header;
  HeadedRecordReader reader(argv[1], &header);
  CHECK(reader.OK());
  // cout << header.DebugString() << endl;

  Show(reader);

  CHECK(reader.OK());

  return 0;
}
