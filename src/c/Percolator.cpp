/**
 * \file Percolator.cpp
 * \brief Object for running percolator
 *****************************************************************************/
#include "Percolator.h"
#include "analyze_psms.h"

using namespace std;


/**
 * \returns a blank Percolator object
 */
Percolator::Percolator() {

}

/**
 * Destructor
 */
Percolator::~Percolator() {
}

/**
 * main method for Percolator
 */
int Percolator::main(int argc, char** argv) {
  analyze_matches_main(PERCOLATOR_COMMAND, argc, argv);
  return 0;
}

/**
 * \returns the command name for Percolator
 */
string Percolator::getName() {
  return "percolator";
}

/**
 * \returns the description for Percolator
 */
string Percolator::getDescription() {
  return 
    "Analyze a collection of PSMs to target and decoy "
    "sequences using the percolator algorithm.";

}

/**
 * \returns the file stem of the application, default getName.
 */
string Percolator::getFileStem() {
  return "percolator";
}

COMMAND_T Percolator::getCommand() {

  return PERCOLATOR_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not.
 */
bool Percolator::needsOutputDirectory() {
  return true;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
