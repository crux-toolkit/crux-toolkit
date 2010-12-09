#include "PrintVersion.h"

#include "version.h"


using namespace std;

PrintVersion::PrintVersion() {

}

PrintVersion::~PrintVersion() {
}


int PrintVersion::main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  printf("Crux version %s\n", VERSION);
  return 0;
}

string PrintVersion::getName() {
  return "version";
}

string PrintVersion::getDescription() {
  return 
    "Print the Crux version number to standard output, "
    "then exit";
}
