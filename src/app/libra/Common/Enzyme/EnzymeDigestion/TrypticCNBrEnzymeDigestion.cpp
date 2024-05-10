#include "TrypticCNBrEnzymeDigestion.h"


/*

Program       : CNBrEnzymeDigestion for PeptideProphet                                                       
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

TrypticCNBrEnzymeDigestion::TrypticCNBrEnzymeDigestion() : EnzymeDigestion("KRM", "", "P", 1, 0) { }

// only count internal M without following modifications
int TrypticCNBrEnzymeDigestion::numMissedCleavages(char* pep) {
  if(strlen(pep) > 4 && pep[1] == '.' && pep[strlen(pep)-2] != '.') {
    cerr << "cannot parse peptide " << pep << endl;
    exit(1);
  }
  int counter = 0;
  int nmc = 0;

  int offset = 0;
  if(strlen(pep) > 2 && pep[1] == '.')
    offset = 2; // have X.PEPTIDE.X format


  Boolean found = False;
  // look for unmodified M's (without symbol afterwords)
  for(int k = offset; k < strlen(pep) - offset; k++) {
    if(counter >= min_dist_ && k >= offset + min_edge_dist_ && k < strlen(pep) - offset - min_edge_dist_) { // position ok
      if(pep[k] == 'M') {
	if(k == strlen(pep) - 1 || (k < strlen(pep) - 1 && pep[k+1] >= 'A' && pep[k+1] <= 'Z')) {
	  nmc++;
	  counter = 0;
	}
      }
      else if(pep[k] == 'K' || pep[k] == 'R') {
	if((k > offset || termCanFollow(pep[k-1])) && (k < strlen(pep) - offset - 1 || termCanPrecede(pep[k+1]))) {
	  nmc++;
	  counter = 0;
	}
      }

    } // ok posn

  } // next pep posn
  return nmc;

}
