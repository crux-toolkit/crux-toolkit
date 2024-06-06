/*
Program       : ASAPRatioPeptideParser
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioProteinRatioParserMain.cpp 7996 2019-12-25 00:16:42Z real_procopio $



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

#include "ASAPRatioProteinRatioParser.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 2) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cerr << " usage: <ASAPRatioProteinRatioParser <protein_xmlfile>" << endl << endl;
    exit(1);
  }

  char *protxml_filename;
  char *testmode=NULL;
  if (argc>2) {
    if (!strncmp(argv[2],REGRESSION_TEST_CMDLINE_ARG,strlen(REGRESSION_TEST_CMDLINE_ARG))) {
        testmode = argv[2];
        protxml_filename = argv[1];
    }
    else {
      testmode = argv[1];
      protxml_filename = argv[2];
    }
  }
  else {
    protxml_filename = argv[1];
  }

  cout << "ASAPRatioProteinRatioParser (" << szTPPVersionInfo << ")" << endl;

  ASAPRatioProteinRatioParser *p = new ASAPRatioProteinRatioParser(protxml_filename,testmode);
  delete p;

  cout << endl << "ASAPRatioProteinRatioParser done." << endl;
  return 0;
}
