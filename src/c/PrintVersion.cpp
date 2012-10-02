/**
 * \file PrintVersion.cpp
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#include "PrintVersion.h"
#include "version.h"
#include "pwiz/Version.hpp"

using namespace std;

/**
 * \returns a blank PrintVersion object
 */
PrintVersion::PrintVersion() {

}

/**
 * Destructor
 */
PrintVersion::~PrintVersion() {
}

/**
 * main method for PrintVersion
 */
int PrintVersion::main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  printf("Crux version %s\n", VERSION);
  
// TODO, get this working (link problem) SJM
//  string proteowizard_version = pwiz::Version::str();
//  printf("Proteowizard version %s\n", proteowizard_version.c_str());

  return 0;
}

/**
 * \returns the command name for PrintVersion
 */
string PrintVersion::getName() {
  return "version";
}

/**
 * \returns the description for PrintVersion
 */
string PrintVersion::getDescription() {
  return 
    "Print the Crux version number to standard output, "
    "then exit.";
}

COMMAND_T PrintVersion::getCommand() {
  return VERSION_COMMAND;
}
