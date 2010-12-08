#include "ComputeQValues.h"
#include "OutputFiles.h"
#include "analyze_psms.h"

using namespace std;

ComputeQValues::ComputeQValues() {

}

ComputeQValues::~ComputeQValues() {
}


int ComputeQValues::main(int argc, char** argv) {
  //TODO : Figure out how to do this.
  analyze_matches_main(QVALUE_COMMAND, argc, argv);
  return 0;
}

string ComputeQValues::getName() {
  return "compute-q-values";
}

string ComputeQValues::getDescription() {
  return 
  "Assign a q-value, which is a statistical confidence "
  "measure that accounts for multiple testing, to each "
  "PSM in a given set";

}

string ComputeQValues::getFileStem() {
  return "qvalues";
}
