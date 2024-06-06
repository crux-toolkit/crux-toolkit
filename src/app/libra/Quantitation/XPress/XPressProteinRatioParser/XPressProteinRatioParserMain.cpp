/*

Program       : XPressProteinRatioParser
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and
                open source code
Date          : 11.27.02
SVN Info      : $Id: XPressProteinRatioParserMain.cpp 7705 2017-12-15 21:43:12Z real_procopio $


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

#include "XPressProteinRatioParser.h"
#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn
#include <fstream>

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths et
  Boolean use_intensities = False;
  char  *testMode = NULL;      // regression test stuff - bpratt Insilicos LLC, Nov 2005

  cout << " " << argv[0] << " (" << szTPPVersionInfo << ")" << endl;

  if (argc < 2) {
    cout << "USAGE:  XPressProteinRatioParser <protxml_file> [options]" << endl;
    cout << "Options:  -i         calculate intensity based protein ratios [default: area-based ratio]" << endl << endl;
    exit(1);
  }

  for (int k = 2; k < argc; k++) {
    if (strlen(argv[k]) > 1 && argv[k][0] == '-') {

      if (argv[k][1] == 'i') {
	use_intensities = True;
	cout << " (calculating intensity-based ratios)" << endl;
      }
      else if (!strncmp (argv[k], REGRESSION_TEST_CMDLINE_ARG, strlen(REGRESSION_TEST_CMDLINE_ARG)))
	testMode = argv[k]; // regression test stuff - bpratt Insilicos LLC, Nov 2005

    }
  }

  //  XPressProteinRatioParser *p = new XPressProteinRatioParser(argv[1],(argc>2)?argv[2]:NULL);
  XPressProteinRatioParser *p = new XPressProteinRatioParser(argv[1], use_intensities, testMode);
  delete p;

  return 0;
}
