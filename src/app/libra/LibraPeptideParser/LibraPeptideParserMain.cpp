#include "LibraPeptideParser.h"
#include "LibraConditionHandler.hpp"
#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn

int main(int argc, char** argv) 
{
  hooks_tpp handler(argc,argv); // set up install paths etc
  char* conditionFileName=NULL;
  char* interactXmlFileName=NULL;

  cout << "LibraPeptideParser (" << szTPPVersionInfo << ")" << endl;

  // Parse arguments
  if ( argc <= 2 ) {
    cout << "usage: LibraPeptideParser <pepXML file> -c<parameter file>"<< endl;
    cout << endl;
    exit(1);
  }

  int fnameArg = 1;
  char *testArg = NULL;

  for ( int argNum = 1 ; argNum < argc ;  argNum++ ) {
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    if (!strncmp(argv[argNum],REGRESSION_TEST_CMDLINE_ARG,
		 strlen(REGRESSION_TEST_CMDLINE_ARG))) {
      // learn or run a test
      testArg = argv[argNum];
      if (1==argNum) {
	fnameArg++; // test arg got in before filename arg
      }
    }
    else {      
      if (argNum == fnameArg) {
	interactXmlFileName = argv[argNum];
      }
      else if( !strncmp( argv[argNum] , "-c", 2 ) ) {
	conditionFileName = argv[argNum] + 2;
      }
    }
  }

  if (conditionFileName==NULL) {
    cout << "Condition File not specified with the -c option!"<< endl;
    cout << "usage: LibraPeptideParser <pepXML file> -c<parameter file>"<< endl;
    cout << endl;
    exit(1);
  }
  else {
    LibraPeptideParser *p = new LibraPeptideParser(interactXmlFileName, conditionFileName, testArg);
    delete p;
  }

  return 0;
}
