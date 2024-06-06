/*
Program       : ASAPRatio
Author        : Xiao-jun Li <xli@systemsbiology.org>
Date          : 09.17.02
SVN Info      : $Id: ASAPRatio_numFns.cpp 7996 2019-12-25 00:16:42Z real_procopio $

Functions for numerical calculation in ASAPRatio.

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "ASAPRatio_numFns.h"
#include "Common/spectStrct.h"

// For a given set of data, dataErrs, and dataWghs, this function identifies 
// any outliers and gets the mean and confidence interval of ratio.
void getDataRatio(double *ratio, double *error, double *h2l_ratio, double *h2l_error, 
		  double confL, 
		  double *data, double *dataErrs, 
		  double *dataWghs, int *dataIndx, int dataSize,
		  int testType)
{

  int counts[4] = {0, 0, 0, 0};
  double *dtRatios, *dtErrors, *dtWeights;
  double *dt_h2l_Ratios,  *dt_h2l_Errors;
  int *dtIndx;
  int pass;
  int count, vldNum;
  double sum;
  double tmpError;
  double acc = 0.01;
  int i, j;


  // check whether there are valid data
  for(i = 0; i < dataSize; ++i) {
    if(data[i] == 0.) 
      ++counts[0];
    else if(data[i] == -1.) 
      ++counts[1];
    else if(data[i] == -2.) 
      ++counts[2];
    else
      ++counts[3];
  }

  // easy output for invalid data set
  if(counts[3] < 1) {
    if(counts[0] > counts[1]) {
      *ratio = 0.;
      *error = 0.;
      *h2l_ratio = 0.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != 0.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    else if(counts[0] < counts[1]) {
      *ratio = -1.;
      *error = 0.;
      *h2l_ratio = -1.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -1.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    else {
      *ratio = -2.;
      *error = 0.;
      *h2l_ratio = -2.;
      *h2l_error = 0.;
      if(testType == 1) {
	for (i = 0; i < dataSize; ++i) 
	  if(data[i] != -2.)
	    dataIndx[i] = 1;
	  else
	    dataIndx[i] = 0;
      }
    }
    return;
  } // if(count[3] < 1) {

  // allocate memory
  dtRatios = (double *) calloc(dataSize, sizeof(double));
  dtErrors = (double *) calloc(dataSize, sizeof(double));
  dt_h2l_Ratios = (double *) calloc(dataSize, sizeof(double));
  dt_h2l_Errors = (double *) calloc(dataSize, sizeof(double));
  dtWeights = (double *) calloc(dataSize, sizeof(double));
  dtIndx = (int *) calloc(dataSize, sizeof(int));


  // identify outliers

  // collect valid data, transform into log(ratio)
  for (i = 0; i < dataSize; ++i) 
    if(data[i] > 0.)
      dataIndx[i] = 0;
    else
      dataIndx[i] = 1;

  count = 0;
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 1) {
	dtRatios[count] = log(data[i]);
	dt_h2l_Ratios[count] = log(1/data[i]); 
	++count;
      }
    }//if(dataIndx[i] == 0) {
  } //for (i = 0; i < dataSize; ++i) {

  // identify any outliers
  if(testType != 0)
    DixonTest(dtRatios, dtIndx, count);

  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) {
      pass = 1;
      for (j = 0; pass == 1 && j < count; ++j) {
	if(dtIndx[j] == 1 
	   && fabs(log(data[i])-dtRatios[j]) < acc*dtRatios[j]){
	  pass = 0;
	}
      }
      if(pass == 0) 
	dataIndx[i] = 1;
    } // if(dataIndx[i] == 0) {
  } //for (i = 0; i < dataSize; ++i) {


  // get ratio and error

  // collect valid date
  count = 0; 
  for (i = 0; i < dataSize; ++i) {
    if(dataIndx[i] == 0) { // not an outlier
      dtRatios[count] = data[i];
      dt_h2l_Ratios[count] = 1 / data[i];
      dtErrors[count] = dataErrs[i];
      dt_h2l_Errors[count] = dataErrs[i] / (data[i] * data[i]);
      dtWeights[count] = dataWghs[i];
      ++count;
    }
  } //for (i = 0; i < dataSize; ++i) {

  // calculate ratio and error
  if(count < 1) { // no valid data
    *ratio = -2.;
    *error = 0.;
    *h2l_ratio = -2.;
    *h2l_error = 0.;
  } //if(count < 1) { // no valid data
  else if(count == 1) { // only one valid data
    *ratio = dtRatios[0];
    *error = dtErrors[0];
    *h2l_ratio =  dt_h2l_Ratios[0];
    *h2l_error =  dt_h2l_Errors[0];
  }//  else if(count == 1) { // only one valid data
  else {
    // transform into log(ratio)
    for (i = 0; i < count; ++i) {
      dtErrors[i] /= dtRatios[i];
      dtRatios[i] = log(dtRatios[i]);
      dt_h2l_Errors[i] /= dt_h2l_Ratios[i];
      dt_h2l_Ratios[i] = log(dt_h2l_Ratios[i]);
     }
    // calculate the light:heavy ratio by weight 
    findMeanAndStdDevWeight(ratio, error, dtRatios, h2l_ratio, h2l_error, dt_h2l_Ratios, dtWeights, count);

    sum = 0.;
    vldNum = 0;
    for (i = 0; i < count; ++i) {
      if(dtErrors[i] > 0.) {
	sum += 1./dtErrors[i]/dtErrors[i];
	++vldNum;
      }
    }
    if(vldNum > 0) {
      tmpError = 1./sqrt(sum);
    } // if(vldNum > 0) {
    else
      tmpError = 0.;

    *error = sqrt((*error)*(*error)+tmpError*tmpError);
    sum = 0.;
    vldNum = 0;
    for (i = 0; i < count; ++i) {
      if(dt_h2l_Errors[i] > 0.) {
	sum += 1./dt_h2l_Errors[i]/dt_h2l_Errors[i];
	++vldNum;
      }
    }
    if(vldNum > 0) {
      tmpError = 1./sqrt(sum);
    } // if(vldNum > 0) {
    else
      tmpError = 0.;

    *h2l_error = sqrt((*h2l_error)*(*h2l_error)+tmpError*tmpError);

    
    // transform back to ratio
    *ratio = exp(*ratio);
    *error *= (*ratio);
    *h2l_ratio = exp(*h2l_ratio);
    *h2l_error *= (*h2l_ratio);
  }//else {
  
  // free memory
  free(dtRatios);
  free(dtErrors);
  free(dtWeights);
  free(dtIndx);

  return;
}


// For a set of data and weight, this function finds the mean and standard deviation.
void findMeanAndStdDevWeight(double *mean, double *error, double *data,
			     double *h2l_mean, double *h2l_error,
			     double *h2l_data, double *weight, int size)
{
  double sum0, sum1, sum2, h2l_sum1, h2l_sum2, sumW;
  double nEff;
  double mnValue, mxValue, mn_inv_Value, mx_inv_Value;
  int count, i;
  
  if (size < 2) {
    *mean = *data;
    *error = 0.;
    *h2l_mean = *h2l_data;
    *h2l_error = 0.;
    return;
  }


  // ensure weight is valid
  count = 0;
  sum0 = 0.;
  mnValue = data[0];
  mxValue = data[0];
  mn_inv_Value = h2l_data[0];
  mx_inv_Value = h2l_data[0];
  for (i = 0; i < size; ++i) {
    if(weight[i] >= 0.) {
      ++count;
      sum0 += weight[i];
    }
    mnValue = mnValue < data[i] ? mnValue : data[i];
    mxValue = mxValue > data[i] ? mxValue : data[i];

    mn_inv_Value = mn_inv_Value < h2l_data[i] ? mn_inv_Value : h2l_data[i];
    mx_inv_Value = mx_inv_Value > h2l_data[i] ? mx_inv_Value : h2l_data[i];
  }

  if(mnValue >= mxValue) {
    *mean = mnValue;
    *error = 0.;
    *h2l_mean = mn_inv_Value;
    *h2l_error = 0.;
    return;
  }

  if(count < size || sum0 == 0.) {   // no all have valid weight
    if(count < 1 || sum0 == 0.) {   // if no data has weight
      for (i = 0; i < size; ++i) {
	weight[i] = 1.;
      }
    } // if(count < 1 || sum0 == 0.) {   // if no data has weight
    else {
      sum0 /= count;
      for (i = 0; i < size; ++i) {
	if(weight[i] < 0.) {
	  weight[i] = sum0;
	}
      }
    }
  } //if(count < size || sum0 == 0.) {   // no all have valid weight

  // get mean and std. dev.
  sum0 = 0.;
  sum1 = 0.;
  sum2 = 0.;
  h2l_sum1 = 0.;
  h2l_sum2 = 0.;
  sumW = 0.;
  for (i = 0; i < size; ++i) {
    sum0 += weight[i];
    sum1 += data[i]*weight[i];
    sum2 += data[i]*data[i]*weight[i];
    h2l_sum1 += h2l_data[i]*weight[i];
    h2l_sum2 += h2l_data[i]*h2l_data[i]*weight[i];
    sumW += weight[i]*weight[i];
  }

  // get mean
  *mean = sum1/sum0;
  *h2l_mean = h2l_sum1/sum0;
  // get std. dev.
  if(sum2*sum0-sum1*sum1 > 0.) {
    nEff = sum0*sum0/sumW;
    if(nEff > 2.) {
      *error = sqrt((sum2*sum0-sum1*sum1)*nEff/(nEff-1.))/sum0;
      *h2l_error = sqrt((h2l_sum2*sum0-h2l_sum1*h2l_sum1)*nEff/(nEff-1.))/sum0;
    }
    else {
      *error = sqrt(2.*(sum2*sum0-sum1*sum1))/sum0;
      *h2l_error = sqrt(2.*(h2l_sum2*sum0-h2l_sum1*h2l_sum1))/sum0;
    }
  }
  else {
    *error = 0.;
    *h2l_error = 0.;
  }
    
  return;
}


// This function uses Dixon's test with alpha = 0.05 to identify any outliers.
void DixonTest(double *data, int *outliers, int size)
{
  double PadeApprx(double x, double *xa, double *ya, int size);

  // cutoff values in Dixon's test: n = 3, ..., 30, INF.
  double ya[29] = {0.941, 0.765, 0.642, 0.560, 0.507, 0.554,
		   0.512, 0.477, 0.576, 0.546, 0.521, 0.546,
		   0.525, 0.507, 0.490, 0.475, 0.462, 0.450,
		   0.440, 0.430, 0.421, 0.413, 0.406, 0.399,
		   0.393, 0.387, 0.381, 0.376, 0.};

  // values of 1/n: 1/3, ..., 1/30, 1/INF.
  double xa[29] = {0.333333, 0.250000, 0.200000, 0.166667, 0.142857, 
		   0.125000, 0.111111, 0.100000, 0.090909, 0.083333, 
		   0.076923, 0.071429, 0.066667, 0.062500, 0.058824, 
		   0.055556, 0.052632, 0.050000, 0.047619, 0.045455, 
		   0.043478, 0.041667, 0.040000, 0.038462, 0.037037, 
		   0.035714, 0.034483, 0.033333, 0.};
  
  int cnstSize = 29;
  int *dataIndx;
  int startIndx, endIndx;
  int count;
  double ratio1, ratio2;
  double cutoff, x;
  int i, j;

  // assume none is an outlier
  for (i = 0; i < size; ++i)
    outliers[i] = 0;
  if (size < 3) // not enough data for checking
    return;
  
  // get dataIndx for ordered data
  dataIndx = (int *) calloc(size, sizeof(int));
  for (i = 0; i < size; ++i)
    dataIndx[i] = i;
  for(i = 0; i < size; ++i) {
    for(j = 0; j < size-i-1; ++j) {
      if(data[dataIndx[j]] > data[dataIndx[j+1]]) {
	count = dataIndx[j];
	dataIndx[j] = dataIndx[j+1];
	dataIndx[j+1] = count;
      }
    }
  } // for(i = 0; i < size; ++i) {
  
  // check for outliers
  count = 0;
  startIndx = 0;
  endIndx = size;
  while(size > 2  // enough data for checking
	&& count != size // look for more when an outlier is identified
	&& data[dataIndx[startIndx]] != data[dataIndx[endIndx-1]]) { 

    // restore size
    count = size;

    // get cutoff
    if (size < 3)
      cutoff = 1.;
    else if(size <= cnstSize+1)
      cutoff = ya[size-3];
    else {
      x = 1./((double) size);
      cutoff = PadeApprx(x, xa, ya, cnstSize);
    }

    // get ratio
    if(size < 8) {
      ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
      ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	/(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx]]);
    }
    else if(size < 11) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+1]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-2]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else if(size < 14) {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-2]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-2]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+1]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+1]]);
      }
      else
	ratio2 = 0.;  
    }
    else {
      if(data[dataIndx[startIndx]] != data[dataIndx[endIndx-3]]) { 
	ratio1 = (data[dataIndx[startIndx+2]]-data[dataIndx[startIndx]])
	  /(data[dataIndx[endIndx-3]]-data[dataIndx[startIndx]]);
      }
      else
	ratio1 = 0.;
      if(data[dataIndx[startIndx+2]] != data[dataIndx[endIndx-1]]) { 
	ratio2 = (data[dataIndx[endIndx-1]]-data[dataIndx[endIndx-3]])
	  /(data[dataIndx[endIndx-1]]-data[dataIndx[startIndx+2]]);
      }
      else
	ratio2 = 0.;  
    }

    // check ratio
    if(ratio1 > ratio2) {
      if(ratio1 > cutoff) { // an outlier
	outliers[dataIndx[startIndx]] = 1;
	--size;
	++startIndx;
      }
    } // if(ratio1 > ratio2) {
    else {
      if(ratio2 > cutoff) { // an outlier
	outliers[dataIndx[endIndx-1]] = 1;
	--size;
	--endIndx;
      }
    } //else {
  } // while(size > 2  // enough data for checking
  

  free(dataIndx);

  return;
}

// This function returns the value of Pade Approximation.
double PadeApprx(double x, double *xa, double *ya, int size)
{
  double y, dy;
  double tiny = 1.e-25;
  int m,i,ns=1;
  double w,t,hh,h,dd,*c,*d;
  double *xb, *yb;
  int n = size;

  // convert into 1 ... n
  xb = (double *) calloc(n+1, sizeof(double));
  yb = (double *) calloc(n+1, sizeof(double));
  for (i = 0; i < size; ++i) {
    xb[i+1] = xa[i];
    yb[i+1] = ya[i];
  }

  // use ratint
  c = (double *) calloc(n+1, sizeof(double));
  d = (double *) calloc(n+1, sizeof(double));

  hh=fabs(x-xb[1]);
  for (i=1;i<=n;i++) {
    h=fabs(x-xb[i]);
    if (h == 0.0) {
      y=yb[i];
      dy=0.0;
      free(c);
      free(d);
      free(xb);
      free(yb);
      return y;
    } else if (h < hh) {
      ns=i;
      hh=h;
    }
    c[i]=yb[i];
    d[i]=yb[i]+tiny;
  }
  y=yb[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      w=c[i+1]-d[i];
      h=xb[i+m]-x;
      t=(xb[i]-x)*d[i]/h;
      dd=t-c[i+1];
      if (dd == 0.0) {
	printf("Error in routine PadeApprx\n");
	free(c);
	free(d);
	free(xb);
	free(yb);
	return y;
      }
      dd=w/dd;
      d[i]=c[i+1]*dd;
      c[i]=t*dd;
    }
    y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  free(c);
  free(d);
  free(xb);
  free(yb);

  return y;
}


// This function converts a ratio and its error into strings for output.
char **ASAPRatio_ratioOutput(double ratio[2], int ratioIndx) 
{
  char **ratioStrings;
  char text[50];
  double numerator, denom;
  int i;

  ratioStrings = (char **) calloc(3, sizeof(char *));
  for (i = 0; i < 3; ++i) {
    ratioStrings[i] = (char *) calloc(10, sizeof(char));
    ratioStrings[i][0] = '\0';
  }

  if(ratio[0] == -2.) { // 0.00 : 0.00
    for (i = 0; i < 3; ++i)
      strcpy(ratioStrings[i], "-1.00");
  }
  else if(ratio[0] == -1.) { // 1.00 : 0.00
    if(ratioIndx == 1) 
      strcpy(ratioStrings[0], "0.00");
    else if(ratioIndx == 2) 
      strcpy(ratioStrings[0], "INF:1");
    else
      strcpy(ratioStrings[0], "9999.");
    for (i = 1; i < 3; ++i)
      strcpy(ratioStrings[i], "0.00");
  }
  else if(ratio[0] == 0.) { // 0.00 : 1.00
    if(ratioIndx == 1) 
      strcpy(ratioStrings[0], "9999.");
    else if(ratioIndx == 2) 
      strcpy(ratioStrings[0], "1:INF");
    else
      strcpy(ratioStrings[0], "0.00");
    for (i = 1; i < 3; ++i)
      strcpy(ratioStrings[i], "0.00");
  }
  else { // normal
    if(ratioIndx == 1) {
      sprintf(ratioStrings[0], "%.2f", 1./ratio[0]);
      sprintf(ratioStrings[1], "%.2f", ratio[1]/ratio[0]/ratio[0]);
    }
    else if(ratioIndx == 2) {
      if(ratio[0] > 1.) {
	sprintf(ratioStrings[0], "%.2f:1", ratio[0]);
	sprintf(ratioStrings[1], "%.2f", ratio[1]);
      }
      else {
	sprintf(ratioStrings[0], "1:%.2f", 1./ratio[0]);
	sprintf(ratioStrings[1], "%.2f", ratio[1]/ratio[0]/ratio[0]);
      }
    }
    else {
      sprintf(ratioStrings[0], "%.2f", ratio[0]);
      sprintf(ratioStrings[1], "%.2f", ratio[1]);
    }

    // want to use only first 2 digits of num and denom for ratio (since xml doesn't store extra digits)
    sprintf(text, "%0.2f", ratio[1]);
    numerator = atof(text);
    sprintf(text, "%0.2f", ratio[0]);
    denom = atof(text);
    if(100.* numerator/denom >= 10)
      sprintf(ratioStrings[2], "%.0f", 100.* numerator/denom);
    else
      sprintf(ratioStrings[2], "%.1f", 100.* numerator/denom);


    //  if(100.*ratio[1]/ratio[0] >= 10)
    //    sprintf(ratioStrings[2], "%.0f", 100.*ratio[1]/ratio[0]);
    //  else
    //    sprintf(ratioStrings[2], "%.1f", 100.*ratio[1]/ratio[0]);
  }

  
  return ratioStrings;
}


// This function frees a matrix.
void freeMtrx(void **mtrx, int size)
{
  int i;
  
  for (i = 0; i < size; ++i)
    free(mtrx[i]);
  free(mtrx);
  
  return;
}
