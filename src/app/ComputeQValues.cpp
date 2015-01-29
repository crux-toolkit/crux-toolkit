/**
 * \file ComputeQValues.cpp
 * \brief Object for calling compute-q-values
 *****************************************************************************/
#include "ComputeQValues.h"
#include "io/OutputFiles.h"
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
  analyze_matches_main(argc, argv);
  return 0;
}

/**
 * \returns the command name for ComputeQValues
 */
string ComputeQValues::getName() {
  return "assign-confidence";
}

/**
 * \returns the description for ComputeQValues
 */
string ComputeQValues::getDescription() {
  return "Assign two types of statistical confidence measures (q-values and "
         "posterior error probabilities) to each PSM in a given set.";
}

/**
 * \returns the filestem for ComputeQValues
 */
string ComputeQValues::getFileStem() {
  return "assign-confidence";
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
