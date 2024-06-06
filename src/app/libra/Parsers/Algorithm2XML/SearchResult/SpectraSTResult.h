#ifndef SPECTRAST_RESULT_H
#define SPECTRAST_RESULT_H

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "SearchResult.h"

/*

Program       : SpectraSTResult for discr_calc of PeptideProphet 
Author        : Henry Lam <hlam@systemsbiology.org>                                                       
Date          : 04.10.06 


Copyright (C) 2006 Andrew Keller

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

Henry Lam
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
hlam@systemsbiology.org

*/

class SpectraSTResult : public SearchResult {

 public:

  SpectraSTResult(Array<Tag*>* tags);
  virtual ~SpectraSTResult(){};

  const char* getName();
  double dot_;
  double delta_;
  double dotbias_;
  double mzdiff_;
  int hitnum_;
  double fval_;
  double libprob_;
  std::string libremark_;

 protected:

 
};











#endif
