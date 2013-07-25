#include "ReadSpectrumRecordsApplication.h"

#include "carp.h"
#include "parameter.h"

ReadSpectrumRecordsApplication::ReadSpectrumRecordsApplication() {
}

ReadSpectrumRecordsApplication::~ReadSpectrumRecordsApplication() {
}

int ReadSpectrumRecordsApplication::main(int argc, char** argv) {

  const char* option_list[] = {
  };
  int num_options = sizeof(option_list) / sizeof(char*);
  const char* arg_list[] = {
    "spectrum records file"
  };
  int num_args = sizeof(arg_list) / sizeof(char*);
  initialize(arg_list, num_args, option_list, num_options, argc, argv);

  carp(CARP_INFO, "Running read-spectrumrecords...");

  string records_file = get_string_parameter_pointer("spectrum records file");

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

string ReadSpectrumRecordsApplication::getName() {
  return "read-spectrumrecords";
}

string ReadSpectrumRecordsApplication::getDescription() {
  return "Runs read-spectrumrecords";
}

bool ReadSpectrumRecordsApplication::needsOutputDirectory() {
  return false;
}

COMMAND_T ReadSpectrumRecordsApplication::getCommand() {
  return READ_SPECTRUMRECORDS_COMMAND;
}

bool ReadSpectrumRecordsApplication::hidden() {
  return true; 
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
