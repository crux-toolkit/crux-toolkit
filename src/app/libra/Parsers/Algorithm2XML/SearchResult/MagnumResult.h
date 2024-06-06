#ifndef Magnum_RESULT_H
#define Magnum_RESULT_H

#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <string>

#include "SearchResult.h"

/*

Program       : MagnumResult for discr_calc of PeptideProphet 
Author        : David Shteynberg <dshteynb%at%systemsbiology.org>                                                       
Date          : 03.28.16


Copyright (C) 2016 David Shteynberg and Andrew Keller

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
#define XL_TYPES 3

class MagnumResult : public SearchResult {

 public:

  MagnumResult(char* szBuf, Boolean preexisting_probs);
  MagnumResult(Array<Tag*>* tags);
  ~MagnumResult();
  void process(char* szBuf, Boolean preex_probs, char* val);
  const char* getName();
  //char* extractDatabase(char* html);

  vector<string>* linkedPeps_;
  vector<vector<string>*>* linkedProts_;
  double magnumScore_;
  double expectScore_;
  double deltaScore_;
  double ppm_error_;

  double sum_rank_;
  
  double lower_score_;
  double higher_score_;

  double alpha_score_;
  double beta_score_;

  double alpha_rank_;
  double beta_rank_;

 protected:

 
};











#endif
