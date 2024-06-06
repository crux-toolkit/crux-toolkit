/*
Program       : Coverage
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02 
SVN Info      : $Id: main.cpp 7712 2017-12-16 01:47:09Z real_procopio $

Program for computing peptide coverage of protein in batch manner from database

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
#include <stdio.h>
#include <stdlib.h>

#include "Coverage.h"


void append(char* orig, char* next, int max_size) {
  int k = strlen(orig);
  int i = 0;
  while(next[i] && k < max_size - 1) {
    if((next[i] >= 'A' && next[i] <= 'Z') || (next[i] >= 'a' && next[i] <= 'z')) {
      orig[k++] = next[i];
    }
    i++;
  }
  // term
  orig[k] = 0;
}

static char *tidy(char *fname) {
#ifdef WINDOWS_NATIVE // single-quote enclosure messes us up...
  for (char *p=fname;p&&*p;) {
    if ('\''==*p) {
      memmove(p,p+1,strlen(p));
    } else {
      p++;
    }
  }
#endif
  return fname;
}


int main(int argc, char** argv) {
  if(argc < 4) {
    cerr << "usage: batchcoverage <database> <inputfile> <outputfile>" << endl;
    exit(1);
  }

  FILE *db = fopen(tidy(argv[1]),"rb");
  if(! db) {
    cerr << "cannot open database " << argv[1] << endl;
    exit(1);
  }

  ifstream fin(tidy(argv[2]));
  if(! fin) {
    cerr << "cannot find coverage file " << argv[2] << endl;
    exit(1);
  }

  ofstream fout(tidy(argv[3]), ios::out);

  const int line_size = 500000;
  char *data =new char[line_size];
  char *next = new char[line_size];
  const int prot_size = 50000;
  char seq[prot_size];
  seq[0] = 0;
  data[0] = 0;
  Coverage* coverage = NULL;
  //Boolean match = False;

  while(fin >> next) {

    if(next[0] == '>') { // new protein
      if(coverage != NULL) {
	fout << coverage->getCoverage() * 100 << endl;
	delete coverage;
	coverage = NULL;
	seq[0] = 0;
      }
      unsigned int nextlen = strlen(next);
      for (int retry=2;retry--;) { // there's an implicit assumption that dbase is alpha sorted
	// if it isn't may need to rewind it
	while(strncmp(next, data, nextlen) && fgets(data, line_size, db)) {
	  strcat(data," ");
	}

	if(!strncmp(next, data, nextlen)) { // match to protein of interest
	  fout << next+1 << "\t";

	  if(fgets(data, line_size,db)) {
	    if(data[0] == '>') {
	      // heave zero length protein
	      fout << -1 << endl;
	    }
	    else {
	      do {
		append(seq, data, prot_size);
	      } while(fgets(data, line_size, db) && (data[0] != '>'));
	      strcat(data," ");
	    }
	  }
	  break; // get out of retry loop
	} // if found next protein
	else if (retry) {
	  data[0] = 0;
	  rewind(db); // rewind
	}
	else {             
	  cerr << "could not find entry for " << next << " in " << argv[1] << endl;
	  exit(1);
	}
      } // end retry loop
    } // if new protein
    else { // a peptide
      if(coverage == NULL) {
	coverage = new Coverage(seq);
	//if(match) {
	//  match = False;
	//  fout << endl << seq << endl;
	//}
      }
      coverage->cover(next);
    }

  } // next coverage data
  // last entry

  if(coverage != NULL) {
    fout << coverage->getCoverage() * 100 << endl;
    delete coverage;
  }
  fclose(db);
  fin.close();
  fout.close();
  delete[] data;
  delete[] next;
  return 0;
}
