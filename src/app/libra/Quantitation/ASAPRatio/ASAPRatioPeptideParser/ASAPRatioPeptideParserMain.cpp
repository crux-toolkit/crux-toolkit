/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: ASAPRatioPeptideParserMain.cpp 7996 2019-12-25 00:16:42Z real_procopio $


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

#include "ASAPRatioPeptideParser.h"
#include "Common/constants.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn
#include <string>

void error(ostream& os) {
  os << "ASAPRatioPeptideParser (" << szTPPVersionInfo <<")" << endl;
  os << "USAGE:  ASAPRatioPeptideParser [xmlfile] [options]" << endl;

  os << "Options:  -l<str>    change labeled residues (default='C')" << endl;
  os << "          -b         heavy labeled peptide elutes before light labeled partner" << endl;
  os << "          -f<num>    areaFlag set to num (ratio display option)" << endl;
  os << "          -r<num>    range around precusor m/z to search for peak (default 0.5)" << endl;
  os << "          -S         static modification quantification (i.e. each run is either all light or all heavy)" << endl;
  os << "          -F         use fixed scan range for light and heavy" << endl;
  os << "          -C         quantitate only the charge state where the CID was made" << endl;
  os << "          -B         return a ratio even if the background is high" << endl;
  os << "          -Z         set all background to zero" << endl;
  os << "          -p<num>    minimum PeptideProphet probability to process (default 0)" << endl;
  os << "          -i<num>    minimum iProphet probability to process (default 0)" << endl;
  os << "          -m<str>    specified label masses (e.g. M74.325Y125.864), only relevant for static modification quantification " << endl;
  os << "          -w         EXPERIMENTAL: use wavelet smoothing from WaveletQuant implmentation" << endl;

  os << endl << "Example1:  ASAPRatioPeptideParser interact"<<get_pepxml_dot_ext()<<" -b" << endl;
  os << endl << "Example2:  ASAPRatioPeptideParser interact"<<get_pepxml_dot_ext()<<" -S -lDE -mD105.23E115.74" << endl;
  os << endl << endl;
  //exit(1);
}

int main(int argc, char** argv) {
  hooks_tpp(argc,argv); // handle install dir issues etc

  InputStruct options;

  char *testMode = NULL; // regression test stuff - bpratt Insilicos LLC, Nov 2005

  // default settings
  options.bUseSameScanRange = TRUE; // -a, -b => (FALSE)
  options.bUseFixedScanRange = FALSE;
  options.bQuantCIDChrgOnly = FALSE;
  options.bQuantHighBackGrnd = FALSE;
  options.bZeroAllBackGrnd = FALSE;
  options.dMassTol = _ASAPRATIO_MZBOUND_;
  options.bXpressLight1 = 0;  // 0 = unused, 1= light (-L), 2=heavy (-H)
  strcpy(options.szXpressResidues, "C"); // -l

  options.bUseWaveletSmoothing = FALSE;

  options.staticQuant = 0;
  options.labelMasses[0] = 0;
  options.dMinPprob = 0;
  options.dMinIprob = 0;

  if(argc < 2) {
    error(cout);
    exit(1);
  }

  // now look for command line
  for(int k = 2; k < argc; k++) {
    if(strlen(argv[k]) > 1 && argv[k][0] == '-')
      if(argv[k][1] == 'b') {
	options.bUseSameScanRange = FALSE;
      }
      else if(argv[k][1] == 'l' && strlen(argv[k]) > 2) {
	strcpy(options.szXpressResidues, argv[k] + 2);
      }
      else if(argv[k][1] == 'f' && strlen(argv[k]) > 2) {
	sscanf(argv[k] + 2, "%d", &options.bXpressLight1);
      }
      else if(argv[k][1] == 'p' && strlen(argv[k]) > 2) {
	sscanf(argv[k] + 2, "%lf", &options.dMinPprob);
      }
      else if(argv[k][1] == 'i' && strlen(argv[k]) > 2) {
	sscanf(argv[k] + 2, "%lf", &options.dMinIprob);
      }
      else if(argv[k][1] == 'r' && strlen(argv[k]) > 2) {
	sscanf(argv[k] + 2, "%lf", &options.dMassTol);
      }
      else if(argv[k][1] == 'S') {
	options.staticQuant = 1;
      }
      else if(argv[k][1] == 'F') {
	options.bUseFixedScanRange = 1;
      }
      else if(argv[k][1] == 'B') {
	options.bQuantHighBackGrnd = 1;
      }
      else if(argv[k][1] == 'Z') {
	options.bZeroAllBackGrnd = 1;
      }
      else if(argv[k][1] == 'C') {
	options.bQuantCIDChrgOnly = 1;
      }
      else if(argv[k][1] == 'w') {
	options.bUseWaveletSmoothing = 1;
      }
      else if(argv[k][1] == 'm' && strlen(argv[k]) > 2) {
	strcpy(options.labelMasses, argv[k] + 2);
      }
      else if (!strncmp(argv[k],REGRESSION_TEST_CMDLINE_ARG,strlen(REGRESSION_TEST_CMDLINE_ARG))) {
	testMode = argv[k]; // regression test stuff - bpratt Insilicos LLC, Nov 2005
      }
      else if(argv[k][1] == 'h') { // display help
	error(cout);

	exit(1);
      }
  }

  cout << "ASAPRatioPeptideParser (" << szTPPVersionInfo <<")" << endl;

  ASAPRatioPeptideParser *p = new ASAPRatioPeptideParser(argv[1], options, testMode, options.dMassTol);
  delete p;

  cout << endl << "ASAPRatioPeptideParser done." << endl;
  return 0;
}
