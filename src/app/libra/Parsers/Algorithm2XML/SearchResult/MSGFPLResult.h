#ifndef MSGFPL_RESULT_H
#define MSGFPL_RESULT_H

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "SearchResult.h"

/*

Program       : MSGFPLResult for discr_calc of PeptideProphet 
Author        : David Shteynberg <dshteynb%at%systemsbiology.org>                                                       
Date          : 09.02.11


Copyright (C) 2007 David Shteynberg and Andrew Keller

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
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

class MSGFPLResult : public SearchResult {

 public:

  MSGFPLResult(char* szBuf, Boolean preexisting_probs);
  MSGFPLResult(Array<Tag*>* tags);
  void process(char* szBuf, Boolean preex_probs, char* val);
  const char* getName();
  //char* extractDatabase(char* html);


  double msgfScore_;
  //  double massError_;
  //double mzSSE_;


 protected:

 
};











#endif
