/**
 * \file PrintProcessedSpectra.cpp
 *
 * DESCRIPTION: Main method for the print-processed-spectra command.
 *              For every spectrum in an ms2 file, process as for
 *              xcorr and print peaks in ms2 format to new file.
 */

#include "PrintProcessedSpectra.h"

#include "print-processed-spectra.h"

using namespace std;

/**
 * \returns a blank PrintProcessedSpectra object
 */
PrintProcessedSpectra::PrintProcessedSpectra() {

}

/**
 * Destructor
 */
PrintProcessedSpectra::~PrintProcessedSpectra() {
}

/**
 * main method for PrintProcessedSpectra
 */
int PrintProcessedSpectra::main(int argc, char** argv) {

  return print_processed_spectra_main(argc, argv);
}

/**
 * \returns the command name for PrintProcessedSpectra
 */
string PrintProcessedSpectra::getName() {
  return "print-processed-spectra";
}

/**
 * \returns the description for PrintProcessedSpectra
 */
string PrintProcessedSpectra::getDescription() {
  return 
    "Write a new ms2 file with all of the same spectra "
    "with only the peaks used for computing xcorr.";
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

