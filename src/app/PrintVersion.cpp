/**
 * \file PrintVersion.cpp
 * \brief Object for printing the crux version number.
 *****************************************************************************/
#include "PrintVersion.h"
#include "crux_version.h"
#include "pwiz/Version.hpp"
#include "Caller.h"
#include "CometSearch/Common.h"
#include "boost/version.hpp"
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

  printf("====================\n");
  printf("Crux version %s\n", CRUX_VERSION);

  printf("====================\n");  
  string proteowizard_version = pwiz::Version::str();
  printf("Proteowizard version %s\n", proteowizard_version.c_str());
  
  printf("====================\n");
  string percolator_version = Caller::greeter();
  printf("%s", percolator_version.c_str());
  
  printf("====================\n");
  printf("Comet version %s\n", comet_version);

  printf("====================\n");
  printf("Boost version %s\n", BOOST_LIB_VERSION);
  printf("====================\n");

  return 0;
}

/**
 * \returns the command name for PrintVersion
 */
string PrintVersion::getName() const {
  return "version";
}

/**
 * \returns the description for PrintVersion
 */
string PrintVersion::getDescription() const {
  return 
    "Print the Crux version number to standard output, "
    "then exit.";
}

COMMAND_T PrintVersion::getCommand() const {
  return VERSION_COMMAND;
}
