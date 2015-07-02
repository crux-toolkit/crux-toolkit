#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <gflags/gflags.h>
#include "spectrum_collection.h"
#include "records.h"
#include "spectrum_preprocess.h"
#include "max_mz.h"

using namespace std;

#define CHECK(x) GOOGLE_CHECK(x)

void DebugPreprocess(RecordReader& reader) {
  MaxMZ::SetGlobalMax(5000);
  for (int specnum = 0; !reader.Done(); ++specnum) {
    pb::Spectrum pb_spectrum;
    reader.Read(&pb_spectrum);
    Spectrum spectrum(pb_spectrum);

    printf("==================== SPECTRUM %d ====================\n", specnum);
    int max_mz;
    ObservedPeakSet observed;
    observed.PreprocessSpectrum(spectrum, 2);
  }
}

void Show(RecordReader& reader) {
  int count = 0;
  while (!reader.Done()) {
    pb::Spectrum pb_spectrum;
    reader.Read(&pb_spectrum);
    Spectrum spectrum(pb_spectrum);
    
    cout << "Spectrum Number: " << spectrum.SpectrumNumber()
    	 << "  Precursor m/z: " << spectrum.PrecursorMZ() << endl;

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

int main(int argc, char* argv[]) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (argc != 2) {
    cerr << "Usage:  " << argv[0] << " SPECTRUM_FILE" << endl;
    return -1;
  }

  RecordReader reader(argv[1]);
  CHECK(reader.OK());

  DebugPreprocess(reader);
  //Show(reader);

  CHECK(reader.OK());

  return 0;
}
