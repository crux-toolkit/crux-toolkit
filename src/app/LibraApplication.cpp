/**
 *\file LibraApplication.cpp 
 *****************************************************************************/
#include "parameter.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "util/MathUtil.h"
#include "LibraApplication.h"
#include <sstream>
#include <iomanip>
#include <ios>

using namespace std;

/**
 * \returns a blank PercolatorApplication object
 */
LibraApplication::LibraApplication() {
}

/**
 * Destructor
 */
LibraApplication::~LibraApplication() {}

/**
 * main method for LibraApplication
 */
int LibraApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("input spectra"));
}

/**
 * \runs Libra application
 */
int LibraApplication::main(const vector<string>& paths) {

  return 0;
}

/**
 * \returns the command name for PercolatorApplication
 */
string LibraApplication::getName() const {
  return "libra";
}

/**
 * \returns the description for PercolatorApplication
 */
string LibraApplication::getDescription() const {
  return
    "[[nohtml:Given a pepXML file,  quantitate the peptides using Libra]]"
    "[[html:<p>Libra is a program that performs quantitation of peptides.</p>]]";
}

/**
 * \returns the command arguments
 */
vector<string> LibraApplication::getArgs() const {
  string arr[] = {
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command options
 */
vector<string> LibraApplication::getOptions() const {
  string arr[] = {
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}

/**
 * \returns the command outputs
 */
vector< pair<string, string> > LibraApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("make-pin.pin",
    "a tab-delimited file containing the input target and decoy PSMs in pin "
    "format. This file can be changed to an absolute path (see --output-file "
    "option)."));
  outputs.push_back(make_pair("make-pin.params.txt",
    "a file containing the name and value of all parameters for the current "
    "operation. Not all parameters in the file may have been used in the "
    "operation. The resulting file can be used with the --parameter-file "
    "option for other crux programs."));
  outputs.push_back(make_pair("make-pin.log.txt",
    "a log file containing a copy of all messages that were printed to "
    "standard error."));
  return outputs;
}

/**
 * \returns whether the application needs the output directory or not. (default false).
 */
bool LibraApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T LibraApplication::getCommand() const {
  return LIBRA_COMMAND;
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
