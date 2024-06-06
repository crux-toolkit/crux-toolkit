/*

Program       : Coverage
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: Coverage.h 7712 2017-12-16 01:47:09Z real_procopio $

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


#ifndef COVERAGE_H
#define COVERAGE_H

#include <string.h>
#include <iostream>
#include <stdlib.h>

#include "sysdepend.h"

class Coverage {

 public:

  Coverage(char* seq);
  ~Coverage();
  float getCoverage();
  void cover(char* pep);
  void cover(char* pep, int offset);

 private:

  char* seq_;
  Boolean* cov_; // paint the coverage

};

using namespace std;

#endif
