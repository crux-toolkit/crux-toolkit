#include "ReadSpectrumRecordsApplication.h"

#include "io/carp.h"
#include "util/Params.h"

ReadSpectrumRecordsApplication::ReadSpectrumRecordsApplication() {
}

ReadSpectrumRecordsApplication::~ReadSpectrumRecordsApplication() {
}

int ReadSpectrumRecordsApplication::main(int argc, char** argv) {
  carp(CARP_INFO, "Running read-spectrumrecords...");

  string records_file = Params::GetString("spectrum records file");

  pb::Header header;
  HeadedRecordReader reader(records_file, &header);
  if (!reader.OK()) {
    carp(CARP_FATAL, "Error reading spectrum records file");
  }

  show(reader);

  if (!reader.OK()) {
    carp(CARP_FATAL, "Error reading spectrum records file");
  }

  return 0;
}

void ReadSpectrumRecordsApplication::show(
  HeadedRecordReader& reader
) {
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

string ReadSpectrumRecordsApplication::getName() const {
  return "read-spectrumrecords";
}

string ReadSpectrumRecordsApplication::getDescription() const {
  return "Runs read-spectrumrecords";
}

vector<string> ReadSpectrumRecordsApplication::getArgs() const {
  string arr[] = {
    "spectrum records file"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

bool ReadSpectrumRecordsApplication::needsOutputDirectory() const {
  return false;
}

COMMAND_T ReadSpectrumRecordsApplication::getCommand() const {
  return READ_SPECTRUMRECORDS_COMMAND;
}

bool ReadSpectrumRecordsApplication::hidden() const {
  return true; 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
