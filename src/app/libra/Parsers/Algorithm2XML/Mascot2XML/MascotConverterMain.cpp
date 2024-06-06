/*

Program       : Mascot2XML converter to pepXML
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


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

#include "MascotConverter.h"
#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Common/util.h"


int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 2) {

     cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << std::endl;
    cerr << " usage: Mascot2XML <file.dat> (-D<database>) (-E<sample_enzyme>) (-html) (-pI) (-notgz) (-nodta) (-desc) (-shortid) (-verbose) (-t(!)<regression test>)" << endl;
    cerr << "     -html : generate html output (default output is pepXML)" << endl;
    cerr << "     -pI : calculate pI for the peptide" << endl;
    cerr << "     -notgz : do not generate (compressed) archive of .dta and .out" << endl;
    cerr << "     -nodta : do not generate .dta files in (compressed) archive; .out files only" << endl;
    cerr << "     -desc : generate protein description in pepXML output" << endl;
    cerr << "     -shortid : use short protein id as per Mascot result (default uses full protein ids in fasta file)" << endl;
    cerr << "     -verbose : tell me more about what is going on" << endl;
		cerr << endl;
    cerr << "   For developers:" << endl;
    cerr << "     use -t!foo.test to create a regression test basis file named foo.test" << endl;
    cerr << "     use -tfoo.test to check against a regression test basis file named foo.test" << endl;
    cerr << "     use -t! or -t to learn/use an automatically named test" << endl;
		cerr << endl;
    cerr << "   Recognised <sample_enyzme> (semi-cleavage specification: semi<sample_enyzme>)" << endl;
		ProteolyticEnzymeFactory enzymeFactory;
		enzymeFactory.showEnzymesCatalog();
		cerr << endl;
    cerr << "   Semi-cleavage can be specified as semi<sample_enyzme>," << endl;
    cerr << "       e.g. semitrypsin, semiargc" << endl;
		cerr << endl;
    cerr << "  ** Please contact spctools-discuss@googlegroups.com" << endl;
    cerr << "     if you wish to add an enzyme missing from the list." << endl;

    exit(1);
  }

  Boolean xml = True; //argc == 4 && ! strcmp(argv[argc-1], "xml");
  Boolean tgz = True;
  Boolean dta = True;
  Boolean description = False;
  Boolean shortProteinId = False;
  Boolean verbose = False;
  char* database = NULL;
  char* sample_enz = NULL;
  Boolean calc_pi = False;
  eTagListFilePurpose testType=NO_TEST;
  std::string testfilename;


  // go through optional args
  for(int k = 2; k < argc; k++) {
    if(strlen(argv[k]) > 1 && argv[k][0] == '-') {
      if(argv[k][1] == 'D') {
        database = new char[strlen(argv[k])-1];
        strcpy(database, argv[k] + 2);
   unCygwinify(database); // no effect in non-cygwin builds
      }
      else if(argv[k][1] == 'E') {
        sample_enz = new char[strlen(argv[k])-1];
        strcpy(sample_enz, argv[k] + 2);
      }
      else if(! strcmp(argv[k]+1, "notgz"))
        tgz = False;
      else if(! strcmp(argv[k]+1, "nodta"))
        dta = False;
      else if(! strcmp(argv[k]+1, "shortid"))
        shortProteinId = True;
      else if(! strcmp(argv[k]+1, "desc"))
        description = True;
      else if(! strcmp(argv[k]+1, "verbose"))
        verbose = True;
      else if(! strcmp(argv[k]+1, "xml"))
        xml = True;
      else if(! strcmp(argv[k]+1, "html"))
        xml = False;
      else if(! strcmp(argv[k]+1, "pI"))
        calc_pi = True;
      else if(!strncmp(argv[k],REGRESSION_TEST_CMDLINE_ARG,strlen(REGRESSION_TEST_CMDLINE_ARG))) { // regression test stuff bpratt Insilicos LLC Dec 2005
        checkRegressionTestArgs(argv[k],testType);
        if (testType!=NO_TEST) {
          char *testFileName=NULL;
          testFileName = constructTagListFilename(argv[1], // input file
             argv[k], // program args
             "Mascot2XML", // program name
             testType); // user info output
          testfilename = testFileName;
          delete[] testFileName;
        }
      }
    } // end if -<arg>
  } // end for k

  MascotConverter *mascot = new MascotConverter(argv[1], xml, database, sample_enz, calc_pi, 
    tgz, dta, description, shortProteinId, verbose,  
    testType, testfilename);

  //
  // regression test stuff bpratt Insilicos LLC Dec 2005
  //
  if (testType != NO_TEST) {
	 const char *progname = argv[0]; // gcc 3.3.2 objects to direct use of argv for some reason
     TagListComparator(progname,testType,mascot->outPepXMLfileName_,testfilename.c_str());
  }

  delete[] database;
  delete[] sample_enz;

  delete mascot;

  return 0;

}
