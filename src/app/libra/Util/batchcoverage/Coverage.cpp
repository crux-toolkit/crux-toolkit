/*

Program       : Coverage
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: Coverage.cpp 7712 2017-12-16 01:47:09Z real_procopio $

Object for computing peptide coverage of protein in batch manner from database

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

#include "Coverage.h"


Coverage::Coverage(char* seq) {
  seq_ = new char[strlen(seq)+1];
  strcpy(seq_, seq);
  seq_[strlen(seq)] = 0;
  cov_ = new Boolean[strlen(seq)];
  for(int k = 0; seq[k]; k++)
    cov_[k] = False;
}

Coverage::~Coverage() {
  if(seq_ != NULL)
    delete [] seq_;
  if(cov_ != NULL)
    delete [] cov_;
}

float Coverage::getCoverage() { 
  int tot = 0;
  for(int k = 0; seq_[k]; k++)
    if(cov_[k])
      tot++;
  //cout << "total matched: " << tot << " out of " << strlen(seq_) << endl;
  return (float)(tot)/strlen(seq_);
}

void Coverage::cover(char* pep) {
  cover(pep, 0);
}

void Coverage::cover(char* pep, int offset) {
  //cout << "pep: " << pep << " offset: " << offset << endl;
  if(offset > (int)(strlen(seq_) - strlen(pep)))
    return;
  char* result = strstr(seq_ + offset, pep);
  if(result != NULL) {
    int diff = strlen(seq_ + offset) - strlen(result);
    if(offset + diff + strlen(pep) > strlen(seq_)) {
      cerr << "error in cover" << endl;
      exit(1);
    }
    //cout << "match from " << offset+diff << " to " << offset + diff + strlen(pep) << endl;
    int limit = offset + diff + strlen(pep);
    for(int k = offset + diff; k < limit; k++)
      cov_[k] = True;
    cover(pep, offset + diff + 1); // continue until no more coverage seen
  }
}
