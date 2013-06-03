#include "TideReadSpectrumRecordsApplication.h"
//#include "tide/read_spectrumrecords.cc"

TideReadSpectrumRecordsApplication::TideReadSpectrumRecordsApplication() {
}

TideReadSpectrumRecordsApplication::~TideReadSpectrumRecordsApplication() {
}

int TideReadSpectrumRecordsApplication::main(int argc, char** argv) {
  return 0;
}

string TideReadSpectrumRecordsApplication::getName() {
  return "tide-read-spectrumrecords";
}

string TideReadSpectrumRecordsApplication::getDescription() {
  return "Runs tide-read-spectrumrecords";
}

bool TideReadSpectrumRecordsApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideReadSpectrumRecordsApplication::getCommand() {
  return TIDE_READ_SPECTRUMRECORDS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
