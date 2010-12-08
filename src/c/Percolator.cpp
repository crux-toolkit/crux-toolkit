#include "Percolator.h"

#include "analyze_psms.h"

using namespace std;

Percolator::Percolator() {

}

Percolator::~Percolator() {
}


int Percolator::main(int argc, char** argv) {
  analyze_matches_main(PERCOLATOR_COMMAND, argc, argv);
  return 0;
}

string Percolator::getName() {
  return "percolator";
}

string Percolator::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the percolator algorithm";

}
