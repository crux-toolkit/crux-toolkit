#include "ReadSpectrumRecordsApplication.h"

#include "io/carp.h"
#include "util/Params.h"
#include "io/SpectrumRecordSpectrumCollection.h"

ReadSpectrumRecordsApplication::ReadSpectrumRecordsApplication() {
}

ReadSpectrumRecordsApplication::~ReadSpectrumRecordsApplication() {
}

int ReadSpectrumRecordsApplication::main(int argc, char** argv) {
  carp(CARP_INFO, "Running read-spectrumrecords...");

  string records_file = Params::GetString("spectrum records file");

  SpectrumRecordSpectrumCollection collection(records_file);
  if (!collection.parse()) {
    carp(CARP_FATAL, "Error parsing spectrum records file");
    return 1;
  }

  for (SpectrumIterator i = collection.begin(); i != collection.end(); i++) {
    cout << "Spectrum Number: " << (*i)->getFirstScan()
         << "  Precursor m/z: " << (*i)->getPrecursorMz() << endl;
    const vector<SpectrumZState>& zStates = (*i)->getZStates();
    cout << "Charge states: " << zStates[0].getCharge();
    for (vector<SpectrumZState>::const_iterator j = zStates.begin() + 1; j != zStates.end(); j++) {
      cout << ", " << j->getCharge();
    }
    cout << endl;
    for (PeakIterator j = (*i)->begin(); j != (*i)->end(); j++) {
      cout << (*j)->getLocation() << " " << (*j)->getIntensity() << endl;
    }
  }

  return 0;
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
