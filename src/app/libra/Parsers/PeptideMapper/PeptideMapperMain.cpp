/*

Program       : PeptideMapper                                                   
Author        : Andrew Keller, Robert Hubley, David Shteynberg and 
                open source code                                                       
Date          : 11.02.2015 

Primary data object holding all mixture distributions for each precursor ion charge

Copyright (C) 2003 Andrew Keller
Copyright (C) 2015 David Shteynberg

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

#include "PeptideMapper.h"

#include "Common/TPPVersion.h" // contains version number, name, revision

#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005

#include <string>

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 3) {
    cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << endl;
    cerr << " usage: PeptideMapper <peptide_file> <protein_database> (<min ntt>(defaults to 0,unless specified in pepXML file with min_number_termini in enzymatic_search_constraint) (DEGEN) (PROT_MW) (PREV_AA_LEN=<length(default=1)>) (NEXT_AA_LEN=<length(default=1)>) (RESTORE_NONEXISTENT_IF_PREFIX=str)" << endl;
#ifdef __LGPL__
    cerr << "Credits: " << "Uses the SPARE Parts library by Bruce W. Watson / Loek Cleophas" << endl;
#endif
    exit(1);
  }
  
  Boolean prot_wt = False;
  int min_ntt = -1;
  Boolean degen = False;
  int n_prev_aas=1;
  int n_next_aas=1;
  std::string testArgs(argv[2]);
  char *restore_prefix = NULL;

  for(int k = 3; k < argc; k++) {
    testArgs += std::string(" ");
    testArgs += std::string(argv[k]);
    if(! strcmp(argv[k], "DEGEN"))
      degen = True;
    else if(! strcmp(argv[k], "PROT_WT"))
      prot_wt = True;
    else if(! strncmp(argv[k], "PREV_AA_LEN=",12)) {
      n_prev_aas = atoi(argv[k]+12);
    } else if(! strncmp(argv[k], "NEXT_AA_LEN=",12)) {
      n_next_aas = atoi(argv[k]+12);
    } else if(! strncmp(argv[k], "RESTORE_NONEXISTENT_IF_PREFIX=",30)) {
      restore_prefix = strdup(argv[k] + 30);
    } else
      min_ntt = atoi(argv[k]);
  }

  PeptideMapper *parser = new PeptideMapper(argv[1], argv[2], min_ntt, n_prev_aas, n_next_aas, degen, prot_wt, testArgs.c_str(), restore_prefix);

  /*
  if(argc == 3)
    new PeptideMapper(argv[1], argv[2], -1, False);
  else if(argc == 4) {
    if(! strcmp(argv[3], "DEGEN"))
      new PeptideMapper(argv[1], argv[2], -1, True);
    else
      new PeptideMapper(argv[1], argv[2], atoi(argv[3]), False);
  }
  else if(argc == 5 && ! strcmp(argv[4], "DEGEN"))
    new PeptideMapper(argv[1], argv[2], atoi(argv[3]), True);
  */
  delete parser;
  return 0;
}
