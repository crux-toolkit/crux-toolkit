// Tandem2xml.cxx
//     Main entry point for converting X!Tandem output to
//     pepXML.

#include <iostream>

using namespace std;

#include "Tandem2xml.h"
#include "TandemResultsParser.h"
#include "Parsers/Algorithm2XML/masscalc.h"
#include "Common/constants.h"

#include "Common/TPPVersion.h" // contains version number, name, revision

int main(int argc, char* argv[]) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if (argc < 2 ) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cerr << "Usage: " << argv[0] << " [USEDESC] [N15] [SCANOFFSETxx] <input-file> [<output-file>] \n";
    cerr << "OPTIONS:\n"
      	 << "\tN15:\tEnable N15 labeling method, (Default: off)\n"
	 << "\tDESCOFF:\tDon't Use Spectrum Descriptions for Naming Spectra in PepXML, (Default: parse scan number from the description)\n"
	 << "\tINDEXOFFxx:\tTPP assumes scans start at 1; older X!Tandem resultshad scan index starts at 0, add xx when converting to pepXML (Default: 0)\n";
    exit(EXIT_FAILURE);
  }

  TandemResultsParser results;

  bool useDescription = true;
  int i = 0;
  int scanOffset = 0;
  bool n15 = false;
  
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "DESCOFF")) {
      useDescription = false;
    }
    else if (!strcmp(argv[i], "N15")) {
      n15 = true;
    }
    else if (!strncmp(argv[i], "INDEXOFF",8)) {
      string arg = string(argv[i] + 8);
      scanOffset = atoi(arg.c_str());
    }
    else {
      break;
    }
  }
  results.setUseN15(n15);
  results.setUseDescription(useDescription);
  results.setIndexOffset(scanOffset);
  results.setFileName(argv[i]);
  if (argc == i+2)
    results.setOutputFile(argv[i+1]);
  
  if (!results.writePepXML())
    exit(EXIT_FAILURE);

  return 0;
}

// Output helpers
string indentChars("   ");
string indentNL("\n");
string& nl()
{return indentNL; }
string& nlIn()
{
  indentNL.append(indentChars);
  return indentNL;
}
string& nlOut()
{
  int lenCur = (int) indentNL.size();
  int lenIndent = (int) indentChars.size();
  if (lenCur > lenIndent)
    indentNL.erase(lenCur - lenIndent);
  return indentNL;
}
