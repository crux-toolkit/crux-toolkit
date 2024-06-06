/*
Program       : PvalueModel
Authors       : Andrew Keller <akeller@systemsbiology.org>
                *Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: PvalueModel.h 7996 2019-12-25 00:16:42Z real_procopio $

Gaussian model used to derive adjusted ratio and stddev as well as pvalue

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

#ifndef PVAL_MOD_H
#define PVAL_MOD_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#include "Common/constants.h"
#include "Quantitation/ASAPRatio/ASAP_structs.h"

class PvalueModel {

 public:
  PvalueModel(double* ratios, int dataNum, const char* pngfile, bool bRegressionTest);
  pValueStrct getParams();

 protected:

  void getDataDistribution(spectStrct *disSpect, double *data, int dataNum);
  void getGaussianFit(double *params, spectStrct &spect, 
		      const char *pngFile, const char *xLbl, const char *yLbl);
  int highestPk(spectStrct &spectrum);
  static double gaussianFit(char **dataStrings, int dataNum,
			    double *params, int pNum);
  
  void get_simplex_minimum(double (*fnctValue)(char **dataStrings, int dataNum,
					       double *params, int pNum),
			   char **dataStrings, int dataNum,
			   double *params, int *paramIndx, int pNum, long *randmSeed);
  
  double simplex_iteration(double (*fnctValue)(char **dataStrings, int dataNum,
					       double *params, int pNum),
			   char **dataStrings, int dataNum,
			   double *params, int *paramIndx, int pNum, 
			   double **simplex, double *dis, 
			   int vldPNum, int imx, double frac);

  double ran2(long *idum);

  pValueStrct pValueData_;

  bool bRegressionTest_; // for regression test use a fixed rand() seed

};

#endif
