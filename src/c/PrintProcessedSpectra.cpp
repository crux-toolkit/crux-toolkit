#include "PrintProcessedSpectra.h"

#include "print-processed-spectra.h"

using namespace std;

PrintProcessedSpectra::PrintProcessedSpectra() {

}

PrintProcessedSpectra::~PrintProcessedSpectra() {
}


int PrintProcessedSpectra::main(int argc, char** argv) {

  return print_processed_spectra_main(argc, argv);
}

string PrintProcessedSpectra::getName() {
  return "print-processed-spectra";
}

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

