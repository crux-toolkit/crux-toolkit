#include "TideResultsApplication.h"
//#include "tide/results.cc"

TideResultsApplication::TideResultsApplication() {
}

TideResultsApplication::~TideResultsApplication() {
}

int TideResultsApplication::main(int argc, char** argv) {
  return 0;
}

string TideResultsApplication::getName() {
  return "tide-results";
}

string TideResultsApplication::getDescription() {
  return "Runs tide-results";
}

bool TideResultsApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideResultsApplication::getCommand() {
  return TIDE_RESULTS_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
