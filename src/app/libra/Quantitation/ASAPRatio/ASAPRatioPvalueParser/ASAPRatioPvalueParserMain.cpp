/*
Program       : main for parser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioPvalueParserMain.cpp 7996 2019-12-25 00:16:42Z real_procopio $

Various XML parsing and overwriting programs

Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "ASAPRatioPvalueParser.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn


int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 2) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cout << " usage: ASAPRatioPvalueParser <xmlfile> [<pngfile>]" << endl;;
    cout << endl << endl;
    exit(1);
  }

  Parser* parser = NULL;
  char *testMode=NULL;
  int filearg = 1;
  int pngarg = 2;
  int nCtorArgs = argc-1;
  for (int i=1;i<argc;i++) {
    if (!strncmp(argv[i],REGRESSION_TEST_CMDLINE_ARG,strlen(REGRESSION_TEST_CMDLINE_ARG))) {
      testMode = argv[i];
      nCtorArgs--;
      if (i==filearg) {
	filearg++;
	pngarg++;
      }
    }
  }
  if(2==nCtorArgs)
    parser = new ASAPRatioPvalueParser(argv[filearg], argv[pngarg], testMode);
  else
    parser = new ASAPRatioPvalueParser(argv[filearg], testMode);
  delete parser;
  return 0;
}
