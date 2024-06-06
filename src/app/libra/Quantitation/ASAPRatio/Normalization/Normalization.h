/*
Program       : Normalization
Authors       : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: Normalization.h 7996 2019-12-25 00:16:42Z real_procopio $

Given Gaussian model parameters, computes adjusted ratio and stddev as well as pvalue

Copyright (C) 2003 Andrew Keller, Xiao-jun Li

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

#ifndef NORMALIZ_H
#define NORMALIZ_H


#include <math.h>
#include <iostream>

#include "Quantitation/ASAPRatio/ASAP_structs.h"

class Normalization {

 public:
  Normalization(pValueStrct params);
  Normalization(double mean, double merr, double stddev);
  double normalize(double* ratio);
  double erfcc(double x);

  double mean_;
  double merr_;
  double stddev_;
}; 


#endif
