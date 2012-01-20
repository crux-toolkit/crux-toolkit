/**
 * \file ComputeQValues.cpp
 * \brief Object for calling compute-q-values
 *****************************************************************************/
#include "ComputeQValues.h"
#include "OutputFiles.h"
#include "analyze_psms.h"

using namespace std;


/**
 * \returns a blank ComputeQValues object
 */
ComputeQValues::ComputeQValues() {

}

/**
 * Destructor
 */
ComputeQValues::~ComputeQValues() {
}

/**
 * main method for ComputeQValues
 */
int ComputeQValues::main(int argc, char** argv) {
  //TODO : Figure out how to do this.
  analyze_matches_main(QVALUE_COMMAND, argc, argv);
  return 0;
}

/**
 * \returns the command name for ComputeQValues
 */
string ComputeQValues::getName() {
  return "compute-q-values";
}

/**
 * \returns the description for ComputeQValues
 */
string ComputeQValues::getDescription() {
  return 
  "Assign a q-value, which is a statistical confidence "
  "measure that accounts for multiple testing, to each "
  "PSM in a given set.";

}

/**
 * \returns the filestem for ComputeQValues
 */
string ComputeQValues::getFileStem() {
  return "qvalues";
}

COMMAND_T ComputeQValues::getCommand() {
  return QVALUE_COMMAND;
}

/**
 * \returns whether the application needs the output directory or not.
 */
bool ComputeQValues::needsOutputDirectory() {
  return true;
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
