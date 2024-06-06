/*
Program       : PvalueModel
Authors       : Andrew Keller <akeller@systemsbiology.org>
                *Xiao-jun Li (xli@systemsbiology.org>
Date          : 11.27.02
SVN info      : $Id: PvalueModel.cpp 7996 2019-12-25 00:16:42Z real_procopio $

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

#include <string>
#include <exception>
#include <unistd.h>

#include "PvalueModel.h"
#include "Common/spectStrct.h"
#include "Common/util.h"

PvalueModel::PvalueModel(double* ratios, int dataNum, const char* pngFile,bool bRegressionTest) {

  double params[6];
  int paramNum = 6;
  spectStrct disSpect;
  double *data;
  int tmpNum;
  int i;

  bRegressionTest_ = bRegressionTest; // for regression tests we fix the rand() seed

  // get log10(uniq_pep_ratio)
  data = (double *) calloc(dataNum, sizeof(double));
  tmpNum = 0;
  for (i = 0; i < dataNum; ++i){
    if(ratios[i] > 0.) {
      data[tmpNum++] = log10(ratios[i]);
    }
  }

  // get distribution
  getDataDistribution(&disSpect, data, tmpNum);
  free(data);

  // get gaussian fit
  getGaussianFit(params, disSpect, pngFile, 
		 "log10(unique L/H peptide ratio)", "unique peptide count");

  // return pValueStrct
  pValueData_.mean = (double)params[1];
  pValueData_.merr = (double)params[4];
  pValueData_.stddev = (double)params[2];
  if(params[0] < 0.)
    pValueData_.merr = -1.;
}

pValueStrct PvalueModel::getParams() { 
  return pValueData_;
}


void PvalueModel::getDataDistribution(spectStrct *disSpect, double *data, int dataNum) {
  double gridSize = 0.02;  // was 0.05 before 2014-05-23
  double lowBound;
  double upBound;
  int gridNum;
  int gridIndx;

  int i;

  lowBound = data[0];
  upBound = data[0];
  for (i = 1; i < dataNum; ++i) {
    if(lowBound > data[i]) lowBound = data[i]; 
    if(upBound < data[i]) upBound = data[i]; 
  }
  lowBound = (floor(lowBound/gridSize+0.5))*gridSize;
  gridNum = (int) ((upBound-lowBound)/gridSize+0.5) + 1;

  disSpect->size = gridNum;
  disSpect->xval = (double *) calloc(disSpect->size, sizeof(double));
  disSpect->yval = (double *) calloc(disSpect->size, sizeof(double));
  for (i = 0; i < disSpect->size; ++i){
    disSpect->xval[i] = lowBound + i*gridSize;
    disSpect->yval[i] = 0.;
  }
  for (i = 0; i < dataNum; ++i) {
    gridIndx = (int) ((data[i]-lowBound)/gridSize+0.5);
    disSpect->yval[gridIndx] += 1.;
  }

  return;
}


void PvalueModel::getGaussianFit(double *params, spectStrct &spect, 
				 const char *pngFile, const char *xLbl, const char *yLbl) {

  // files
  FILE *fout;
  std::string tmpFile1(pngFile);
  std::string tmpFile2(pngFile);
  std::string tmpFile3(pngFile);

  int paraNum = 6;
  int paramIndx[6] = {1, 1, 1, 0, 0, 0};
  double fitDis = 0.01;

  char **dataStrings;
  int dataNum = spect.size;
  spectStrct smSpect, fitData, origData;
  double gridSize;
  int dtIndx, pkIndx;
  double x, dx;
  double f, df[3], ddf[3][3];
  double **matrx, col[3];
  int indx[3];
  double error, zValue;
  long randmSeed = 0;
  time_t currTime;
  int psIndx = 0;
  double tmpParams[6], initParams[6];
  int ite = 10;
  double mnValue=-1, tmpValue;
 
  int i, j, k;

  // check gridSize
  if(spect.size < 2){
    printf("Not enough data for Gaussian fitting. \n");
    params[0] = -1.;
    return;
  }
  else
    gridSize = spect.xval[1] - spect.xval[0];

  // get randMSeed
  (void) time(&currTime);
  srand48((long) (bRegressionTest_?1234:currTime));
  randmSeed = -lrand48();

  origData = spect;

  // get initial guess
  if(pngFile != NULL) {
    tmpFile1 += "tmp_raw_data.dat.XXXXXX"; 
    // do this in designated tmp dir, if any
    replace_path_with_webserver_tmp(tmpFile1);

    if((fout = FILE_mkstemp(tmpFile1)) == NULL) {  
      printf("Cannot write to \"%s\".\n", tmpFile1.c_str());
      params[0] = -1.;
      return;
    }
    for (i = 0; i < spect.size; ++i)
      fprintf(fout, "%f\t%f\n", spect.xval[i], spect.yval[i]);
    fclose(fout);
  }

  // smooth
  smSpect = spect;
  smSpect.smoothSpectFlx(5*gridSize, 4);

  // get fine-tuned curve
  fitData.size = (int)(smSpect.xval[smSpect.size-1]/fitDis) - (int)(smSpect.xval[0]/fitDis) + 10;
  fitData.xval = (double *) calloc(fitData.size, sizeof(double));
  fitData.yval = (double *) calloc(fitData.size, sizeof(double));
  x = ((int)(smSpect.xval[0]/fitDis)+1)*fitDis; 
  fitData.size = 0;
  dtIndx = 0;
  while(x < smSpect.xval[smSpect.size-1]) {
    fitData.xval[fitData.size] = x;
    while(dtIndx < smSpect.size && x > smSpect.xval[dtIndx])
      ++dtIndx;
    fitData.yval[fitData.size] = smSpect.yval[dtIndx-1] 
      + (x-smSpect.xval[dtIndx-1])*(smSpect.yval[dtIndx]-smSpect.yval[dtIndx-1])
      /(smSpect.xval[dtIndx]-smSpect.xval[dtIndx-1]);
    ++fitData.size;
    x += fitDis;
  }
  smSpect.clear();
  fitData.smoothSpectFlx(5*gridSize, 4);

  // get initial guess
  if((pkIndx = highestPk(fitData)) != -1){
    params[0] = fitData.yval[pkIndx];
    params[1] = fitData.xval[pkIndx];
    dtIndx = pkIndx;
    while(dtIndx >= 0 && fitData.yval[dtIndx] > 0.5*fitData.yval[pkIndx])
      --dtIndx;
    if(dtIndx < 0){
      //printf("Cannot fit with gaussian curve 1.\n");
      params[0] = -1.;
      verified_unlink(tmpFile1);
      return;
    }
    else {
      params[2] = -fitData.xval[dtIndx];
    }

    dtIndx = pkIndx;
    while(dtIndx < fitData.size && fitData.yval[dtIndx] > 0.5*fitData.yval[pkIndx])
      ++dtIndx;
    if(dtIndx >= fitData.size) {
      //printf("Cannot fit with gaussian curve 2.\n");
      params[0] = -1.;
      verified_unlink(tmpFile1);
      return;
    }
    else {
      params[2] += fitData.xval[dtIndx];
    }
    params[2] /= sqrt(8.*log(2.)); 
  }
  else {
    //printf("Cannot fit with gaussian curve 3.\n");
    params[0] = -1.;
    verified_unlink(tmpFile1);
    return;
  }
  fitData.clear();

  // gaussian fit
  // remove far away data from fitting
  dataNum = 0;
  for(i = 0; i < spect.size; ++i) {
    if(fabs(spect.xval[i]-params[1]) < 4.*params[2]) {
      spect.xval[dataNum] = spect.xval[i];
      spect.yval[dataNum] = spect.yval[i];
      ++dataNum;
    }
  }
  spect.size = dataNum;

  // collect dataStrings
  dataStrings = (char **) calloc(spect.size, sizeof(char *));
  for (i = 0; i < spect.size; ++i){
    dataStrings[i] = (char *) calloc(1000, sizeof(char));
    sprintf(dataStrings[i], "%f\t%f\n", spect.xval[i], spect.yval[i]);
  }
  for (j = 0; j < paraNum; ++j)
    initParams[j] = params[j];
  for (i = 0; i < ite; ++i){
    for (j = 0; j < paraNum; ++j)
      tmpParams[j] = initParams[j];
    get_simplex_minimum(gaussianFit, dataStrings, dataNum, 
			tmpParams, paramIndx, paraNum, &randmSeed);
    tmpValue = gaussianFit(dataStrings, dataNum,
			   tmpParams, paraNum);
    if(i == 0 || tmpValue < mnValue){
      mnValue = tmpValue;
      for (j = 0; j < paraNum; ++j)
	params[j] = tmpParams[j];
    }
  }
  for(i = 0; i < dataNum; ++i)
    free(dataStrings[i]);
  free(dataStrings);

  // estimate errors in params
  // estimate error in data
  error = 0.;
  for (i = 0; i < spect.size; ++i){
    zValue = (spect.xval[i]-params[1])/params[2];
    error += (spect.yval[i]-params[0]*exp(-0.5*zValue*zValue))
      *(spect.yval[i]-params[0]*exp(-0.5*zValue*zValue));
  }
  if(spect.size > 3.)
    error = sqrt(error/(spect.size-3.));
  else
    error = 1.;

  // get matrix
  matrx = (double **) calloc(3, sizeof(double *));
  for (i = 0; i < 3; ++i){
    matrx[i] = (double *) calloc(3, sizeof(double));
    for (j = 0; j < 3; ++j)
      matrx[i][j] = 0.;
    col[i] = 0.;
  }
  for (i = 0; i < spect.size; ++i) {
    zValue = (spect.xval[i]-params[1])/params[2];
    f = params[0]*exp(-0.5*zValue*zValue);
    df[0] = f/params[0];
    df[1] = f*zValue/params[2];
    df[2] = f*zValue*zValue/params[2];
    ddf[0][0] = 0.;
    ddf[0][1] = df[1]/params[0];
    ddf[1][0] = ddf[0][1];
    ddf[0][2] = df[2]/params[0];
    ddf[2][0] = ddf[0][2];
    ddf[1][1] = f*(zValue*zValue-1.)/params[2]/params[2];
    ddf[1][2] = f*(zValue*zValue-2.)*zValue/params[2]/params[2];
    ddf[2][1] = ddf[1][2];
    ddf[2][2] = f*(zValue*zValue-3.)*zValue*zValue/params[2]/params[2];

    for (j = 0; j < 3; ++j){
      for(k = j; k < 3; ++k){
	matrx[j][k] += 2.*spect.yval[i]*spect.yval[i]*df[j]*df[k]/error/error
	  - 2.*spect.yval[i]*spect.yval[i]*zValue*ddf[j][k]/error;
	matrx[k][j] = matrx[j][k];
      }
    }
  }

  // reverse matrix
  myLUDcmp(matrx, 3, indx, &x);

  // get fit error on params[1]
  for (i = 0; i < 3; ++i)
    params[3+i] = 0.;
  
  for (i = 0; i < spect.size; ++i) {
    zValue = (spect.xval[i]-params[1])/params[2];
    f = params[0]*exp(-0.5*zValue*zValue);
    df[0] = f/params[0];
    df[1] = f*zValue/params[2];
    df[2] = f*zValue*zValue/params[2];
    for (j = 0; j < 3; ++j) 
      col[j] = 2.*spect.yval[i]*spect.yval[i]*df[j]/error/error; 
    myLUBksb(matrx, 3, indx, col);
    for (j = 0; j < 3; ++j)
      params[3+j] += col[j]*col[j]*error*error;
  }
  for (i = 0; i < 3; ++i)
    params[3+i] = sqrt(params[3+i]);
  for (i = 0; i < 3; ++i)
    free(matrx[i]);
  free(matrx);

  // generate graph
  if(pngFile != NULL) {
    tmpFile2 += "tmp_gaussian_fit.dat.XXXXXX"; 
    // do this in designated tmp dir, if any
    replace_path_with_webserver_tmp(tmpFile2);

    if((fout = FILE_mkstemp(tmpFile2)) == NULL) {  
      printf("Cannot write to \"%s\".\n", tmpFile2.c_str());
      verified_unlink(tmpFile1);
      return;
    }
    x = params[1] - 5.*params[2];
    dx = 0.02*params[2];
    while(x <= params[1] + 5.*params[2]){
      fprintf(fout, "%f\t%f\n", x, 
	      params[0]*exp(-0.5*(x-params[1])*(x-params[1])/params[2]/params[2]));
      x += dx;
    }
    fclose(fout);

    if(strstr(pngFile, ".ps") != NULL)
      psIndx = 1;
    tmpFile3 += ".gnuplot.XXXXXX";
    // do this in designated tmp dir, if any
    replace_path_with_webserver_tmp(tmpFile3);

    if((fout = FILE_mkstemp(tmpFile3)) == NULL) {  
      printf("Cannot open file \"%s\".\n",tmpFile3.c_str());
      verified_unlink(tmpFile1);
      verified_unlink(tmpFile2);
      return;
    }
    fprintf(fout, "set grid\n");
    fprintf(fout, "set xlabel \"%s\"\n", xLbl);
    fprintf(fout, "set ylabel \"%s\"\n", yLbl);
    fprintf(fout, "set label \"Amplitude: %.3f+-%.3f\"at graph 0.05, 0.85 left\n", params[0], params[3]);
    fprintf(fout, "set label \"Mean Ratio: %.3f+-%.3f\"at graph 0.05, 0.8 left\n", 
	    pow(10., params[1]), params[4]*log(10.)*pow(10., params[1]));
    if(psIndx != 1) {
      fprintf(fout, "set label \"Ratio Std Dev: %.3f+-%.3f\"at graph 0.05, 0.75 left\n", params[2], params[5]);
      fprintf(fout, "set terminal png\n");
    }
    else {
      fprintf(fout, "set label \"Ratio Std Dev: %.3f+-%.3f\"at graph 0.05, 0.75 left\n", params[2], params[5]);
      fprintf(fout, "set t post enhanced 20\n");
    }
    fprintf(fout, "set o \'%s\'\n", pngFile);
    fprintf(fout, "pl \'%s\'t \"data\"w l lw 3, \'%s\'t \"fitting\"w l lt 3 lw 3\n", tmpFile1.c_str(), tmpFile2.c_str());
    fprintf(fout, "unset output\n");
    fprintf(fout, "quit\n");
    fflush(fout);
    fclose(fout);
    std::string gnuCommand(GNUPLOT_BINARY);
    gnuCommand += " ";
    gnuCommand += tmpFile3;

    try {
      verified_unlink(pngFile);
      verified_system(gnuCommand.c_str()); // system() with verbose error check
    }
    catch (const std::exception& e) {
    }

    //Sleep Max 30 secs
    FILE* test;
    int count = 0;
    while (count++ < 30 && (test = fopen(pngFile,"r"))==NULL) {
      usleep(1000);
    }
    if (count >= 30) {
      fprintf(stderr,"WARNING: Max Timeout reached waiting for gnuplot ... ");
    }
    else if ( test!=NULL ) {
      fclose(test);
    }

    //printf("tmpFile1=%s, tmpFile2=%s, pngFile=%s\n", tmpFile1, tmpFile2, pngFile);

    verified_unlink(tmpFile1);
    verified_unlink(tmpFile2);
    verified_unlink(tmpFile3);

    // write model data for use in protXML
    char* tmpFile4 = new char[strlen(pngFile)+12];
    strcpy(tmpFile4, pngFile);
    strcat(tmpFile4, ".model.data");

    if ((fout = fopen(tmpFile4,"w"))!=NULL) {
      for (i = 0; i < origData.size; ++i) {

	fprintf(fout, "<point logratio=\"%f\" obs_distr=\"%f\" model_distr=\"%f\"/>\n",
		origData.xval[i],
		origData.yval[i],
		params[0]*exp(-0.5*(origData.xval[i]-params[1])*(origData.xval[i]-params[1])/params[2]/params[2])
		);
      }
      fclose(fout);
    }

    delete[] tmpFile4;
  }

  return;
}

int PvalueModel::highestPk(spectStrct &spectrum) {
  int pkIndx=-1;
  double intensity = 0.;
  int i;

  for (i = 0; i < spectrum.size; ++i) {
    if (spectrum.yval[i] > intensity) {
      pkIndx = i;
      intensity = spectrum.yval[i];
    }
  }

  return pkIndx;
}


double PvalueModel::gaussianFit(char **dataStrings, int dataNum,
				double *params, int pNum) {
  double diff, zValue;
  double *xVal, *yVal;
  int i;

  xVal = (double *) calloc(dataNum, sizeof(double));
  yVal = (double *) calloc(dataNum, sizeof(double));
  for (i = 0; i < dataNum; ++i) {
    sscanf(dataStrings[i], "%lf %lf", xVal+i, yVal+i);
  }

  // get diff;
  diff = 0.;
  for (i = 0; i < dataNum; ++i) {
    zValue = (xVal[i]-params[1])/params[2];
    diff += (yVal[i]-params[0]*exp(-0.5*zValue*zValue))
      *(yVal[i]-params[0]*exp(-0.5*zValue*zValue))*yVal[i];
  }
  free(xVal);
  free(yVal);

  return diff;
}


void PvalueModel::get_simplex_minimum(double (*fnctValue)(char **dataStrings, int dataNum,
							  double *params, int pNum),
				      char **dataStrings, int dataNum,
				      double *params, int *paramIndx, int pNum, long *randmSeed) {
  double acc = 1.e-9;
  double tol = 1.e-9; // tolerance on flatness of function
  double dp = 0.2;
  int vldPNum;
  double **simplex; // (ndim+1, ndim+1) matrix for simplex: simplex[][ndim] 
	           // is for the corresponding functional value 
  double *dis;
  double frac;
  int imx, imn, inmx;
  double ymx, ymn, ynmx;
  double yTry;
  int pCount1, pCount2;
  double sum, rtol, distance, tiny = 1.e-10;
  int iteration;
  int i, j, k;

  // get # of valid params
  vldPNum = 0;
  for (i = 0; i < pNum; ++i) {
    if(paramIndx[i] == 1)
      ++vldPNum;
  }

  // check whether there are valid parameters
  if(vldPNum == 0) {
    printf("No valid paramters.\n");
    return;
  }

  // get simplex
  simplex = (double **) calloc(vldPNum+1, sizeof(double *));
  simplex[0] = (double *) calloc(vldPNum+1, sizeof(double));
  pCount1 = 0;
  for (i = 0; i < pNum; ++i) {
    if(paramIndx[i] == 1) {
      simplex[0][pCount1] = params[i];
      ++pCount1;
    }
  }
  yTry = (*fnctValue)(dataStrings, dataNum, params, pNum);
  simplex[0][vldPNum] = yTry;

  pCount1 = 1;
  for (i = 0; i < pNum; ++i) {
    if(paramIndx[i] == 1) {
      simplex[pCount1] = (double *) calloc(vldPNum+1, sizeof(double));
      frac = 1. + dp*(ran2(randmSeed)-0.5);
      params[i] *= frac;
      pCount2 = 0;
      for (j = 0; j < pNum; ++j) {
	if(paramIndx[j] == 1) {
	  simplex[pCount1][pCount2] = params[j];
	  ++pCount2;
	}
      }
      yTry = (*fnctValue)(dataStrings, dataNum, params, pNum);
      simplex[pCount1][vldPNum] = yTry;
      params[i] /= frac;
      ++pCount1;
    }
  }

  //find solution iteratively
  dis = (double *) calloc(vldPNum, sizeof(double));
  iteration = 0;
  while(iteration < 10000) {
    // find ymx, ynmx, ymn
    if(simplex[0][vldPNum] > simplex[1][vldPNum]){
      imx = 0;
      ymx = simplex[0][vldPNum];
      inmx = 1;
      ynmx = simplex[1][vldPNum];
    }
    else {
      imx = 1;
      ymx = simplex[1][vldPNum];
      inmx = 0;
      ynmx = simplex[0][vldPNum];
    }
    imn = 0;
    ymn = simplex[0][vldPNum];

    for (i = 0; i <= vldPNum; ++i) {
      yTry = simplex[i][vldPNum];
      if(yTry > ymx) {
	ynmx = ymx;
	inmx = imx;
	ymx = yTry;
	imx = i;
      }
      else if(yTry > ynmx && i != imx) {
	ynmx = yTry;
	inmx = i;
      }
      else if(yTry < ymn) {
	ymn = yTry;
	imn = i;
      }
    }

    // get params minimize fnctValue
    pCount1 = 0;
    for(i = 0; i < pNum; ++i) {
      if(paramIndx[i] == 1) {
	params[i] = simplex[imn][pCount1];
	++pCount1;
      }
    }

    // get distance bewteen vertex
    distance = 0.;
    for (i = 0; i <= vldPNum; ++i){
      for (j = i+1; j <= vldPNum; ++j){
	sum = 0.;
	for (k = 0; k < vldPNum; ++k){
	  sum += (simplex[i][k]-simplex[j][k])*(simplex[i][k]-simplex[j][k]);
	}
	sum = sqrt(sum);
	if(sum > distance)
	  distance = sum;
      }
    }

    // functional flatnes
    rtol = 2.*fabs(ymx-ymn)/(fabs(ymx)+fabs(ymn)+tiny);

    // check for convergence
    if(distance < acc && rtol < tol) {
      iteration = -1;
      //      printf("\n");
      break;
    }

    // get displacement
    for (i = 0; i < vldPNum; ++i) {
      sum = 0.;
      for (j = 0; j <= vldPNum; ++j) {
	sum += simplex[j][i] - simplex[imx][i];
      }
      dis[i] = sum;
    }

    // reflection
    yTry = simplex_iteration(fnctValue, dataStrings, dataNum, 
			     params, paramIndx, pNum, 
			     simplex, dis, vldPNum, imx, 1.); 

    if(yTry <= ymn) { // extrapolation
      yTry = simplex_iteration(fnctValue, dataStrings, dataNum, 
			       params, paramIndx, pNum, 
			       simplex, dis, vldPNum, imx, -1.); 
    }
    else if(yTry >= ynmx) {
      ymx = simplex[imx][vldPNum];
      // contraction
      yTry = simplex_iteration(fnctValue, dataStrings, dataNum, 
			       params, paramIndx, pNum, 
			       simplex, dis, vldPNum, imx, 0.25); 
      if(yTry >= ymx){ //multiple contraction
	for(i = 0; i <= vldPNum; ++i) {
	  if (i == imn) 
	    continue;
	  for(j = 0; j < vldPNum; ++j) {
	    dis[j] = 0.5*(simplex[i][j]+simplex[imn][j]);
	    simplex[i][j] = dis[j];
	  }

	  pCount1 = 0;
	  for(j = 0; j < pNum; ++j) {
	    if(paramIndx[j] == 1) {
	      params[j] = dis[pCount1];
	      ++pCount1;
	    }
	  }
	  yTry = (*fnctValue)(dataStrings, dataNum, params, pNum);
	  simplex[i][vldPNum] = yTry;
	}
      }
    }

    ++iteration;
  }
  free(dis);
  for(i = 0; i <= vldPNum; ++i)
    free(simplex[i]);
  free(simplex);

  return;
}


double PvalueModel::simplex_iteration(double (*fnctValue)(char **dataStrings, int dataNum,
							  double *params, int pNum),
				      char **dataStrings, int dataNum,
				      double *params, int *paramIndx, int pNum, 
				      double **simplex, double *dis, 
				      int vldPNum, int imx, double frac) {

  double *ptry;
  double ytry;
  int pCount;
  int i;

  // get ptry
  ptry = (double *) calloc(vldPNum+1, sizeof(double));
  for(i = 0; i < vldPNum; ++i) {
    ptry[i] = simplex[imx][i] + frac*dis[i];
  }

  // convert to params
  pCount = 0;
  for(i = 0; i < pNum; ++i) {
    if(paramIndx[i] == 1) {
      params[i] = ptry[pCount];
      ++pCount;
    }
  }
  ytry = (*fnctValue)(dataStrings, dataNum, params, pNum);
  ptry[vldPNum] = ytry;

  // replace if needed
  if(ytry < simplex[imx][vldPNum]) {
    for(i = 0; i <= vldPNum; ++i) {
      simplex[imx][i] = ptry[i];
    }
    for(i = 0; i < vldPNum; ++i) 
      dis[i] *= -1.;
  }

  free(ptry);

  return ytry;
}


double PvalueModel::ran2(long *idum){
  long IM1 = 2147483563;
  long IM2 = 2147483399;
  double AM = 1.0/IM1;
  long IMM1 = IM1-1;
  int IA1 = 40014;
  int IA2 = 40692;
  int IQ1 = 53668;
  int IQ2 = 52774;
  int IR1 = 12211;
  int IR2 = 3791;
  int NTAB = 32;
  long NDIV = (1+IMM1/NTAB);
  double EPS = 1.2e-7;
  double RNMX = (1.0-EPS);

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[32];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;

}
