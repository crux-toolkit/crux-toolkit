#include "TideIndexApplication.h"
//#include "tide/index.cc"

TideIndexApplication::TideIndexApplication() {
}

TideIndexApplication::~TideIndexApplication() {
}

int TideIndexApplication::main(int argc, char** argv) {
  return 0;
}

string TideIndexApplication::getName() {
  return "tide-index";
}

string TideIndexApplication::getDescription() {
  return "Runs tide-index";
}

bool TideIndexApplication::needsOutputDirectory() {
  return true;
}

COMMAND_T TideIndexApplication::getCommand() {
  return TIDE_INDEX_COMMAND;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
