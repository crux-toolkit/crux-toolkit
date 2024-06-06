/*
Program       : ASAPRatio
Author        : Xiao-jun Li <xli@systemsbiology.org>
Date          : 09.17.02
SVN Info      : $Id: ASAPRatio_numFns.h 7996 2019-12-25 00:16:42Z real_procopio $

Header file of ASAPRatio_numFns.c

Copyright (C) 2002 Xiao-jun Li

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

Xiao-jun Li
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
xli@systemsbiology.org
*/

#ifndef _ASAPRATIO_NUMFNS_H_
#define _ASAPRATIO_NUMFNS_H_


//
// Constants
//
#define _ASAPRATIO_CONFL_ 0.6826 // confidence level at 1 sigma

#include "Common/spectStrct.h"


//
// Functions
//

// For a given set of data, dataErrs, and dataWghs, this function identifies 
// any outliers and gets the mean and confidence interval of ratio.
void getDataRatio(double *ratio, double *error, double *h2l_ratio, double *h2l_error, 
		  double confL, double *data, double *dataErrs, 
		  double *dataWghs, int *dataIndx, int dataSize,
		  int testType);

// This function converts a ratio and its error into strings for output.
char **ASAPRatio_ratioOutput(double ratio[2], int ratioIndx);

// This function frees a matrix.
void freeMtrx(void **mtrx, int size);

// This function uses Dixon's test with alpha = 0.05 to identify any outliers.
void DixonTest(double *data, int *outliers, int size);

// For a set of data and weight, this function finds the mean and standard deviation.
void findMeanAndStdDevWeight(double *mean, double *error, double *data, double *h2l_mean,
			     double *h2l_error, double *h2l_data, double *weight, int size);

#endif /* _ASAPRATIO_NUMFNS_H_ */
