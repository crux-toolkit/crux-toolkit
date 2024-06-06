#ifndef COMRESULT_H
#define COMRESULT_H

#include<string.h>
#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>


#include "SearchResult.h"
#include "Parsers/Parser/Tag.h"
#include "Common/constants.h"
#include "Common/sysdepend.h"


/*

Program       : CometResult for discr_calc of PeptideProphet 
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


class CometResult : public SearchResult {

 public:

  CometResult(char* szBuf, Boolean preexistingprobs, int format);
  CometResult(Array<Tag*>* tags);
  CometResult();
  void process(char* szBuf, Boolean preexisting_probs, int format, int* formats, int numformats, char* valid);
  Boolean parseUnknown(const char* szBuf, Boolean preexisting_probs, int* formats, int numformats);
  virtual Boolean parse(const char* szBuf, Boolean preexisting_probs, int format);
  virtual char* extractDatabase(char* html);
  const char* getName();
  double xcorr_;
  double delta_;
  double deltastar_;
  int rank_;
  double sp_score_;
  double expect_;

 protected:

};

#endif
