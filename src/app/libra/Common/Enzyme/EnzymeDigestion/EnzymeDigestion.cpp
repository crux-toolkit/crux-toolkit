#include "EnzymeDigestion.h"

/*

Program       : EnzymeDigestion for PeptideProphet                                                       
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

EnzymeDigestion::EnzymeDigestion(char* sites, char* term_not_following, char* term_not_preceding, 
				 int min_edge_dist, int min_dist) {
  recognition_sites_ = strCopy(sites);
  term_not_following_ = strCopy(term_not_following);
  term_not_preceding_ = strCopy(term_not_preceding);
  min_edge_dist_ = min_edge_dist;
  min_dist_ = min_dist;
}

EnzymeDigestion::~EnzymeDigestion() {
  if(recognition_sites_ != NULL)
    delete recognition_sites_;
  if(term_not_following_ != NULL)
    delete term_not_following_;
  if(term_not_preceding_ != NULL)
    delete term_not_preceding_;
}


int EnzymeDigestion::numMissedCleavages(char* peptide) {
  char* pep = strip(peptide, 1);
  if(strlen(pep) > 4 && pep[1] == '.' && pep[strlen(pep)-2] != '.') {
    cerr << "cannot parse peptide " << pep << endl;
    exit(1);
  }
  int MIN_PEP_LENGTH = 9;
  if(strlen(pep) < MIN_PEP_LENGTH)
    return -1;

  int counter = 0;
  int nmc = 0;
  //nmc = 0.0;
  Boolean found = False;
  int offset = 0;
  if(strlen(peptide) > 2 && peptide[1] == '.')
    offset = 2; // have X.PEPTIDE.X format


  for(int k = offset; k < strlen(pep) - offset; k++) {
    if(counter >= min_dist_ && k >= min_edge_dist_ + offset && k < strlen(pep) - offset - min_edge_dist_) { // position ok
      if((k > offset || termCanFollow(pep[k-1])) && (k < strlen(pep) - offset - 1 || termCanPrecede(pep[k+1])) &&
	 isCompatibleTerminus(pep[k])) {
	nmc++;
	counter = 0;
      }
    } // ok posn

  } // next pep posn
  if(pep != NULL)
    delete pep;
  return nmc;
  return (int)((double)(nmc) * 100 / (strlen(pep) - min_edge_dist_));

}

/////////////////////////////////////////////////////////////
// This was implemented directly 9.30.04
int EnzymeDigestion::numCompatibleTermini(char prev, char* pep, char foll) {
  if(prev == '?' || foll == '?')
    return 2;

  int nct = 0;
  if(prev == '-' || prev == '1' || prev == 'M') // M in case at start of protein (after M cleaved off)
    nct++;
  else if(strlen(pep) > 0 && termCanPrecede(pep[0]))
    for(int s = 0; s < strlen(recognition_sites_); s++) 
      if(prev == recognition_sites_[s]) {
	nct++;
	s = strlen(recognition_sites_); // done
      }
  if(foll == '1' || foll == '-')
    nct++;
  else if((strlen(pep) > 1 && termCanFollow(pep[strlen(pep)-2])) || 
	  (strlen(pep) == 1 && termCanFollow(prev)))
    for(int s = 0; s < strlen(recognition_sites_); s++) 
      if(pep[strlen(pep)-1] == recognition_sites_[s]) {
	nct++;
	s = strlen(recognition_sites_); // done
      }
      
  // cout << prev << " " << pep << " " << foll << ": " << nct << endl;
  return nct;
  /*



  char* worker = new char[strlen(pep)+5];
  worker[0] = prev;
  worker[1] = '.';
  worker[2] = 0;
  strcat(worker, pep);
  worker[strlen(pep)+2] = '.';
  worker[strlen(pep)+3] = foll;
  worker[strlen(pep)+4] = 0;
  int result = numCompatibleTermini(worker);
  if(worker != NULL)
    delete worker;
  return result;
  */
}



// pep of form: P.XXXXX.F
int EnzymeDigestion::numCompatibleTermini(char* pep) {
  if(strlen(pep) < 4 || pep[1] != '.' || pep[strlen(pep)-2] != '.') {
    cerr << "cannot parse peptide " << pep << endl;
    exit(1);
  }
  int nct = 0;
  if(pep[0] == '-' || pep[0] == '1' || pep[0] == 'M') // M in case at start of protein (after M cleaved off)
    nct++;
  else if(strlen(pep) > 2 && termCanPrecede(pep[2]))
    for(int s = 0; s < strlen(recognition_sites_); s++) 
      if(pep[0] == recognition_sites_[s]) {
	nct++;
	s = strlen(recognition_sites_); // done
      }
  if(pep[strlen(pep)-1] == '1' || pep[strlen(pep)-1] == '-')
    nct++;
  else if(strlen(pep) > 3 && termCanFollow(pep[strlen(pep)-4]))
    for(int s = 0; s < strlen(recognition_sites_); s++) 
      if(pep[strlen(pep)-3] == recognition_sites_[s]) {
	nct++;
	s = strlen(recognition_sites_); // done
      }

  return nct;
}


char* EnzymeDigestion::strip(char* pep, Boolean remove_mods) {
  int start = 0;
  int stop = strlen(pep)-1;
  char* output = NULL;
  if(strlen(pep) > 4 && pep[1] == '.')
    start = 2;
  if(strlen(pep) > 4 && pep[strlen(pep)-2] == '.')
    stop = strlen(pep)-3;

  if(! remove_mods) {
    output = new char[stop - start + 2];
    strncpy(output, pep+start, stop - start + 1);
    output[stop - start + 1] = 0;
    return output;
  }

  int length = 0;
  for(int k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      length++;
  output = new char[length+1];
  length = 0;
  for(int k = start; k <= stop; k++) 
    if(pep[k] >= 'A' && pep[k] <= 'Z')
      output[length++] = pep[k];
  output[length] = 0;
  return output;
}


Boolean EnzymeDigestion::isCompatibleTerminus(char c) {
  for(int k = 0; k < strlen(recognition_sites_); k++)
    if(recognition_sites_[k] == c)
      return True;
  return False;
}

Boolean EnzymeDigestion::termCanPrecede(char c) {
  for(int k = 0; k < strlen(term_not_preceding_); k++)
    if(term_not_preceding_[k] == c)
      return False;
  return True;
}

Boolean EnzymeDigestion::termCanFollow(char c) {
  for(int k = 0; k < strlen(term_not_following_); k++)
    if(term_not_following_[k] == c)
      return False;
  return True;
}

char* EnzymeDigestion::strCopy(const char* orig) {
  char* output = new char[strlen(orig)+1];
  strcpy(output, orig);
  output[strlen(orig)] = 0;
  return output;
}
