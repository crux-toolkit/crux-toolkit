/*
Program       : Normalization
Author        : Andrew Keller <akeller@systemsbiology.org>
                Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: Normalization.cpp 7996 2019-12-25 00:16:42Z real_procopio $

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

#include "Normalization.h"


Normalization::Normalization(pValueStrct params) {
  mean_ = params.mean;
  merr_ = params.merr;
  stddev_ = params.stddev;
}

Normalization::Normalization(double mean, double merr, double stddev) {
  mean_ = mean;
  merr_ = merr;
  stddev_ = stddev;
}  

double Normalization::normalize(double* ratio) {

  //cout << "mean: " << mean_ << " merr: " << merr_ << endl;
  //cout << "raio: " << ratio[0] << "+/-" << ratio[1] << endl;

  if(merr_ < 0.0)
    return -1.0;
 
  if(ratio[0] > 999.0)
    return 0.0; // zero pvalue
  else if(ratio[0] < 0.0)
    return -1.0; // illegal value, ignore

  double nRatio[2], tmpRatio[2];
  double relErr, z, pVal;
  double min_pvalue = 0.00000001; // 10-8

  // get normalization factor
  tmpRatio[0] = pow(10., (double) mean_);
  tmpRatio[1] = log(10.)*tmpRatio[0]*merr_;

  if(ratio[0] > 0.){
    nRatio[0] = ratio[0]/tmpRatio[0];      
    nRatio[1] = nRatio[0]
      *sqrt(ratio[1]*ratio[1]/ratio[0]/ratio[0]
	    +tmpRatio[1]*tmpRatio[1]/tmpRatio[0]/tmpRatio[0]);
    
    relErr = ratio[1]/ratio[0]/log(10.);
    z = fabs((log10(ratio[0])-mean_)
	     /sqrt(2.*(relErr*relErr+stddev_*stddev_
		       +merr_*merr_)));
    pVal = erfcc(z);
  }
  else if(ratio[0] == 0.){
    nRatio[0] = 0.;
    nRatio[1] = 0.;
    pVal = 0.;
  }
  else if(ratio[0] == -1.) {
    nRatio[0] = -1.;
    nRatio[1] = 0.;
    pVal = 0.;
  }
  else {
    nRatio[0] = -2.;
    nRatio[1] = 0.;
    pVal = -1.;
  }

  ratio[0] = nRatio[0];
  ratio[1] = nRatio[1];

  if(ratio[0] > 999.)
    ratio[0] = 999.; // max possible value

  // put in the min value here
  if(pVal < min_pvalue)
    return 0.0;
  return pVal;
}

double Normalization::erfcc(double x) {
  double t,z,ans;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	    t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}
