/*

Program       : Quantic
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>
Date          : 04.20.2018
SVN Info      : $Id$

Copyright (C) 2018 David Shteynberg

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA3

*/

#include "QuanticParser.h"
#include "Common/util.h"
#include "Common/TPPVersion.h"
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn


int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc
  
  //TODO: Add input error handling
  if(argc < 2) {
    //    cerr <<  argv[0] << " (" << szTPPVersionInfo << ")" << endl;
    cerr << "USAGE: " << argv[0] << " <OPTIONS> <input_file"<<get_pepxml_dot_ext()<<"> [<output_file>]"<< endl 
	 << "\nOPTIONS\n"
	 << "\tTOPNPEAKS=<number>\t Use specified number of peaks for antic score, set to 0 to use all (default=6)\n"
         << "\tMZTOL=<number>\t Use specified +/- MS2 mz tolerance on site specific ions (default=0.1 m/z)\n"
	 << "\tMINPROB=<number>\t Use specified minimum probability to evaluate peptides (default=0)\n"
      	 << "\tANNOTATE\t Report a string of all matched ions used in summation (default=off)\n"
      	 << "\tCOUPLED\t For quantitation use also peptide coupled neutral losses (default=off, use only modified fragment neutral losses)\n"

	 << "\tDIAMODE\t Correct for Isotopes in DIA Mode (default=off, do not correct)\n"

	 << "\t<amino_acid_token1><mass_shift1>,<amino_acids_token2><mass_shift2>...\tSpecify molecular formulas for modifications (e.g. \"K[300]C7H12N2O3,M[147]C5H9NO2S,n[29]N2H\" )\n"
	 << "\tMAXTHREADS=<number>\t Use specified number of threads for processing (default=1)\n"
       << "\t<amino acids, n, or c>:<mass_shift>:<neut_loss1>:...:<neut_lossN>\tSpecify modifications with neutral losses for quantitation (default=off )\n"
         << "\nVERSION: " << szTPPVersionInfo 
	 << endl ;
    exit(1);
  }

  // regression test stuff- bpratt Insilicos LLC
  int testarg = 0;
  int max_threads = 1;
  unsigned int verbose = 0;
  eTagListFilePurpose testType=NO_TEST;
  char *testArgArg=NULL;
  char *testFileName=NULL;
  string* catFile = NULL;
  string* arg = NULL;

  string optString = "";
  double minProb =-100;
  double rtMinProb = 0.9;
  double mztol = -1;
  double ppmtol = 1;


  string ptmMolForms = "";
  
  double massdiff_offset = 0.;

  bool massdiff_mode = false;

  bool dia_mode = false;

  bool lability = false;
  bool direct = false;
  bool autodirect = false;

  bool pep_coupled = false;
  bool annotate = false;
  string modString = "";

  unsigned int top_peaks = 6;
  int em = 1;
  bool update = true;
  bool mpx = true;
  bool keepold = false;

  if (!strncmp(argv[1],REGRESSION_TEST_CMDLINE_ARG,strlen(REGRESSION_TEST_CMDLINE_ARG))) {
     checkRegressionTestArgs(testArgArg = strdup(argv[testarg=1]),testType);
     if (testType!=NO_TEST) {  
        testFileName = constructTagListFilename(argv[1+testarg], // input data file
           testArgArg, // args to the program
           "QuanticParser",  // program name
           testType); // user info output
    }
  }

  for (int argidx = testarg+1; argidx < argc-1; argidx++) {
    string arg = argv[argidx];
    
    if (!strncmp(argv[argidx], "MAXTHREADS=", 11)) {
      arg = string(argv[argidx] + 11);
      max_threads = atoi(arg.c_str());
    }
    else if (!strncmp(argv[argidx], "TOPNPEAKS=", 10)) {
      arg = string(argv[argidx] + 10);
      top_peaks = atoi(arg.c_str());
    }
    else if (!strcmp(argv[argidx], "NOUPDATE")) {
      update = false;
    }
    else if (!strcmp(argv[argidx], "KEEPOLD")) {
      keepold = true;
    }
    else if (!strcmp(argv[argidx], "DIAMODE")) {
      dia_mode = true;
    }
    else if (!strcmp(argv[argidx], "ANNOTATE")) {
      annotate = true;
    }
    else if (!strcmp(argv[argidx], "COUPLED")) {
      pep_coupled = true;
    }
    else if (strstr(argv[argidx], "]")!=NULL) {
      ptmMolForms = arg;
    }
    else if (!strncmp(argv[argidx], "MZTOL=", 6)) {
      arg = string(argv[argidx] + 6);
      mztol = atof(arg.c_str());
    }
    else if (!strncmp(argv[argidx], "MINPROB=", 8)) {
      arg = string(argv[argidx] + 8);
      minProb = atof(arg.c_str());
    }
    else if (strstr(argv[argidx], ":")!=NULL && strstr(argv[argidx], "/")==NULL && strstr(argv[argidx], "\\")==NULL) {
      modString = arg;
    }



    optString += string(argv[argidx]) + ' ' ;
  }

  optString += string(argv[argc-1]);


  QuanticParser* p;
  

  p = new QuanticParser(mztol, ppmtol, verbose);

  if (!modString.empty()) {
    p->setModString(modString);
  }

  if (!ptmMolForms.empty()) {
    p->setMolForms(ptmMolForms);

  }
    
  p->setAnnotate(annotate);
  p->setUpdate(update);
  p->setMinProb(minProb);
  p->setPeptideCoupled(pep_coupled);
  p->setDiaMode(dia_mode);
  
  string in_file = "";
  string out_file = "";
  for (int argidx = testarg+1; argidx<argc; argidx++) {
    if ( strcmp(argv[argidx], "COUPLED") && strcmp(argv[argidx], "DIAMODE") && strcmp(argv[argidx], "ANNOTATE") && strcmp(argv[argidx], "NOUPDATE") && strncmp(argv[argidx], "MZTOL=",6) && strncmp(argv[argidx], "PPMTOL=",7) && strcmp(argv[argidx], "KEEPOLD") && strncmp(argv[argidx], "MINPROB=",8) &&  strncmp(argv[argidx], "MAXTHREADS=",11)  &&  strncmp(argv[argidx], "TOPNPEAKS=",10) && 	  strstr(argv[argidx], "]")==NULL &&
	 ( strstr(argv[argidx], ":")==NULL || 
	   strstr(argv[argidx], "/")!=NULL || 
	   strstr(argv[argidx], "\\")!=NULL ) ) {
      
      if (in_file.empty()) {
	in_file = argv[argidx];
      }
      else {
	out_file = argv[argidx];
      }
    }
    
  }

  if (!out_file.empty()) {
    p->setOutFile(out_file);
  }
  else {
    out_file = in_file;
  }
  p->run(in_file.c_str(), optString.c_str(), max_threads, top_peaks);
  
 
  const char *outfilename = out_file.c_str();
  // regression test?
  if (testType!=NO_TEST) {
    TagListComparator("QuanticParser",testType,outfilename,testFileName);
    free(testArgArg);
    delete[] testFileName;
  }
}
