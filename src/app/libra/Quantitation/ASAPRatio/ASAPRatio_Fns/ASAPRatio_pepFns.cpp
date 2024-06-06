/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatio_pepFns.cpp 8386 2021-02-09 22:48:16Z dshteyn $


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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <sys/param.h>
#endif
#include <sys/stat.h> 
#ifdef __CYGWIN__
#include <sys/dir.h>
#elif !defined(WIN32)
//#include <sys/dir.h>
#endif

#include "Common/constants.h"
#include "Common/spectStrct.h"
#include "ASAPRatio_numFns.h"
#include "ASAPRatio_pepFns.h"

#include "Parsers/mzParser/mzParser.h"
#include "Common/util.h"
#include "Common/sysdepend.h"
#include <string>

using namespace mzParser;

static int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
		     double timeWd, char *xmlFile);
static int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
		     double timeWd, char *xmlFile, double mzBound, int smoothItrs, double smoothRTwindow);
static int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, 
			  double timeWd, char *xmlFile, double mzBound, int smoothItrs, double smoothRTwindow, bool wavelet);
static void getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, 
		      const spectStrct &fitSpectrum, int scan, int scanRange);
static void getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, 
		      const spectStrct &fitSpectrum, int scan, int scanRange, int areaType, Boolean quantHighBG, Boolean zeroBG);
static int getLCSpectAreaAndTime(double *area, double *peakTime, 
			    const spectStrct &rawSpectrum, 
			    const spectStrct &fitSpectrum, 
			    int pkScan, int *valleyScan, double background);
static int getLCSpectAreaAndTime(double *area, double *peakTime, 
			    const spectStrct &rawSpectrum, 
			    const spectStrct &fitSpectrum, 
			    int pkScan, int *valleyScan, double background, int areaType, Boolean quantHighBG);
static void plotPepPeak(const char *pngFileBase, const pepDataStrct &data, 
		 const lcSpectStrct &lcSpect, int qState); 
static void jsPlotPepPeak(const pepDataStrct &data, const lcSpectStrct &lcSpect, int qState); 
static double getLCSpectBackground(const spectStrct &rawSpectrum, 
			      const spectStrct &fitSpectrum, 
			      int *ranges, int *valleyScan);
static double getLCSpectBackground(const spectStrct &rawSpectrum, 
			      const spectStrct &fitSpectrum, 
			    int *ranges, int *valleyScan, Boolean zeroBG);
static int spectPeak(const spectStrct &spectrum, int position, int direction);
static int findNextValley(const spectStrct &spectrum, int position, int direction);
static void getLCSpectPeakAndValleys(int *peakScan, int *valleyScan, 
			      const spectStrct &spectrum, 
			      int scan, double background);


/*
  This function evaluates a peptide ratio from peptide sequence.
  It returns NULL for invalid data.

  peptide specific variables: char *pepSeq, long scan (end scan), int chrg, char *xmlFile 
  global variables: 
  int eltn, int areaFlag, residueStrct *modAAs, int modAANum, pairStrct *prtnAAs, int prtnAANum

*/
pepDataStrct *evalPepDataStrct(char *pepSeq, long scan, int chrg, char *xmlFile,
			       int eltn, int areaFlag, residueStrct *modAAs, int modAANum,
			       pairStrct *prtnAAs, int prtnAANum)
{
  void getPairSequences(char *sequence, char *prtnSequence, 
			int *cidIndx, double *msLight, double *msHeavy, 
			residueStrct *modAAs, int modAANum,
			pairStrct *prtnAAs, int prtnAANum);
  void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		       char *pngFileBase, int cgiIndx);

  pepDataStrct *data;
  char prtnSeq[1000];
  double msLight, msHeavy;
  int cidIndx;

  // initialize
  data = (pepDataStrct *)calloc(1, sizeof(pepDataStrct));
  data->indx = -1;
  
  // light and heavy sequences
  getPairSequences(pepSeq, prtnSeq, &cidIndx, &msLight, &msHeavy, 
		   modAAs, modAANum, prtnAAs, prtnAANum);
  
  // get peptide data
  if(msLight <= 0. || msHeavy <= 0. 
     || msHeavy <= msLight) {
    free(data);
    return NULL;
  }

  // initialize data
  data->scan = scan;
  data->chrg = chrg;
  data->cidIndx = cidIndx;
  data->msLight = msLight;
  data->msHeavy = msHeavy;
  data->eltn = eltn;
  data->areaFlag = areaFlag;

  // calculate pepDataStrct
  getPepDataStrct(data, xmlFile, NULL, 0);

  if(data->indx < 0){
    free(data);
    return NULL;
  }
  else
    return data;
}


// This function returns a light sequence and its heavy partner.
void getPairSequences(char *sequence, char *prtnSequence, 
		      int *cidIndx, double *msLight, double *msHeavy, 
		      residueStrct *modAAs, int modAANum,
		      pairStrct *prtnAAs, int prtnAANum)
{
  double getSeqMonoMZ(char *sequence, int qState,
		      residueStrct *modifiedAA, int modAANum);

  char tmpSequence[1000];
  int seqIndx, prtnIndx;
  int lngth = (int)strlen(sequence);
  char end[3];
  int flag;
  double msTh, msPt;
  int i;


  // get partner sequence

  // check for N-end
  seqIndx = 0;
  prtnIndx = 0;
  if(sequence[0] == '+') { // heavy N-end
    strcpy(end, "<+");
    ++seqIndx;
  }
  else{ // not heavy N-end
    strcpy(end, "<");
  }
  flag = 0; // check whether has a partner
  for (i = 0; i < prtnAANum; ++i) {
    if(strncmp(end, prtnAAs[i].prtnA, strlen(prtnAAs[i].prtnA)) == 0
       || strncmp(end, prtnAAs[i].prtnB, strlen(prtnAAs[i].prtnB)) == 0) {
      flag = 1;
      break;
    }
  }
  if(flag == 1) { // has a partner
    if(sequence[0] != '+') { // not heavy N-end
      prtnSequence[prtnIndx++] = '+';
    }
  } //   if(flag == 1) { // has a partner
  else { // no partner
    if(sequence[0] == '+') { // heavy N-end
      prtnSequence[prtnIndx++] = '+';
    }
  } //   else { // no partner

  // check sequence one by one
  while(seqIndx < lngth) {
    flag = -1; // assume having no partner
    // check on prtnA
    for (i = 0; i < prtnAANum && flag == -1; ++i) {
      if(strncmp(sequence+seqIndx, prtnAAs[i].prtnA, 
		 strlen(prtnAAs[i].prtnA)) == 0) {
	flag = i;
      }
    } // for (i = 0; i < prtnAANum; ++i) {
    if (flag != -1) { // has a partner
      strcpy(prtnSequence+prtnIndx, prtnAAs[flag].prtnB);
      seqIndx += (int)strlen(prtnAAs[flag].prtnA);
      prtnIndx += (int)strlen(prtnAAs[flag].prtnB);
      continue;
    }
    // check on prtnB
    for (i = 0; i < prtnAANum && flag == -1; ++i) {
      if(strncmp(sequence+seqIndx, prtnAAs[i].prtnB, 
		 strlen(prtnAAs[i].prtnB)) == 0) {
	flag = i;
      }
    } // for (i = 0; i < prtnAANum; ++i) {
    if (flag == -1) { // has no partner
      prtnSequence[prtnIndx++] = sequence[seqIndx++];
    }
    else { // has a partner
      strcpy(prtnSequence+prtnIndx, prtnAAs[flag].prtnA);
      seqIndx += (int)strlen(prtnAAs[flag].prtnB);
      prtnIndx += (int)strlen(prtnAAs[flag].prtnA);
    }
  } // while(seqIndx < lngth) {

  // check for C-end
  if(sequence[lngth-1] == '-') { // heavy C-end
    strcpy(end, ">-");
    --prtnIndx;
  }
  else{ // not heavy C-end
    strcpy(end, ">");
  }
  flag = 0; // check whether has a partner
  for (i = 0; i < prtnAANum; ++i) {
    if(strncmp(end, prtnAAs[i].prtnA, strlen(prtnAAs[i].prtnA)) == 0
       || strncmp(end, prtnAAs[i].prtnB, strlen(prtnAAs[i].prtnB)) == 0) {
      flag = 1;
      break;
    }
  }
  if(flag == 1) { // has a partner
    if(sequence[lngth-1] != '-') { // not heavy C-end
      prtnSequence[prtnIndx++] = '-';
    }
  } //   if(flag == 1) { // has a partner
  else { // no partner
    if(sequence[lngth-1] == '-') { // heavy C-end
      prtnSequence[prtnIndx++] = '-';
    }
  } //   else { // no partner

  prtnSequence[prtnIndx] = '\0';


  // decide sequence

  // theoretical mass
  msTh = getSeqMonoMZ(sequence, 0, modAAs, modAANum);
  msPt = getSeqMonoMZ(prtnSequence, 0, modAAs, modAANum);
  if(msTh <= 0. || msPt <= 0. || msTh == msPt){
    *msLight = -1.;
    *msHeavy = -1.;
    if(strcmp(sequence, prtnSequence) > 0) {
      *cidIndx = 1;
      strcpy(tmpSequence, sequence);
      strcpy(sequence, prtnSequence);
      strcpy(prtnSequence, tmpSequence);
    }
    else
      *cidIndx = 0;    
  }
  else if(msTh < msPt){
    *cidIndx = 0;
    *msLight = msTh;
    *msHeavy = msPt;
  }
  else if(msTh > msPt){
    *cidIndx = 1;
    *msLight = msPt;
    *msHeavy = msTh;
    strcpy(tmpSequence, sequence);
    strcpy(sequence, prtnSequence);
    strcpy(prtnSequence, tmpSequence);
  }

  return;
}


// For a given sequence "char * sequence" of charge "int qState", this function gives its (mono-isotopic mass)/charge.
double getSeqMonoMZ(char *sequence, int qState, residueStrct *modifiedAA, int modAANum)
{
  double mass;
  int length = (int)strlen(sequence);
  int count = 0; 
  int flag;
  char end[3];
  char *tmpSeq;
  int i;

  // only positive or neutral charge states allowed
  if (qState < 0) {
    printf("Cannot handle negative charge state.\n");
    return -1.;
  }

  // check input
  getRidOfSpace(sequence);
  cnvtUpper(sequence);
  
  // mass of "qState*H"
  mass = qState*nativeAA[0].ms;


  // mass of N-end
  if(sequence[0] != '+') { // not heavy N-end
    strcpy(end, "<");
  }
  else{ // heavy N-end
    ++count;
    strcpy(end, "<+");
  }

  flag = 0;
  // check for modified Amino Acids first
  for (i = 0; i < modAANum && flag == 0; ++i) {
    if(strncmp(end, modifiedAA[i].rp, modifiedAA[i].sz) == 0) {
      flag = 1;
      mass += modifiedAA[i].ms;
    }
  }
  // check for native Amino Acids
  for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
    if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
      flag = 1;
      mass += nativeAA[i].ms;
    }
  }
  // if no match
  if(flag == 0) {
    printf("N-end unspecified.\n");
    return -1.;
  }


  // mass of C-end
  if(sequence[length-1] != '-') { // not heavy C-end
    strcpy(end, ">");
  }
  else{
    strcpy(end, ">-");
    --length;
  }

  flag = 0;
  // check for modified Amino Acids first
  for (i = 0; i < modAANum && flag == 0; ++i) {
    if(strncmp(end, modifiedAA[i].rp, modifiedAA[i].sz) == 0) {
      flag = 1;
      mass += modifiedAA[i].ms;
    }
  }
  // check for native Amino Acids
  for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
    if(strncmp(end, nativeAA[i].rp, nativeAA[i].sz) == 0) {
      flag = 1;
      mass += nativeAA[i].ms;
    }
  }
  // if no match
  if(flag == 0) {
    printf("C-end unspecified.\n");
    return -1.;
  }

  // add sequence mass
  while(count < length) {
    flag = 0;

    // check for modified Amino Acids first
    for (i = 0; i < modAANum && flag == 0; ++i) {
      if(strncmp(sequence+count, modifiedAA[i].rp, modifiedAA[i].sz) == 0) {
	flag = 1;
	mass += modifiedAA[i].ms;
	count += modifiedAA[i].sz;
      }
    }
    // check for native Amino Acids
    for (i = 0; i < NATIVEAANUM && flag == 0; ++i) {
      if(strncmp(sequence+count, nativeAA[i].rp, nativeAA[i].sz) == 0) {
	flag = 1;
	mass += nativeAA[i].ms;
	count += nativeAA[i].sz;
      }
    }
    // if no match
    if(flag == 0) {
      tmpSeq = (char *) calloc(count+1, sizeof(count));
      strncpy(tmpSeq, sequence, count);
      tmpSeq[count] = '\0';
      printf("Unknown amino acid at \"%s'%c'%s\". \n", 
	     tmpSeq, sequence[count], sequence+count+1);
      free(tmpSeq);
      return -1.;
    }
  } // while(count < length) {
  
  if(qState > 0)
    return mass/qState;
  else
    return mass;
}


// This function evaluates the light:heavy and light:heavy ratio of a peptide.
// For cgi display: cgiIndx = 1 and pngFileBase != NULL.
// For js/flot    : cgiIndx = 2
void getPepDataStrct(pepDataStrct *data, char *xmlFile, char *pngFileBase, int cgiIndx) {
  void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		       char *pngFileBase, int cgiIndx, 
		       Boolean quantHighBG, Boolean zeroBG, 
		       double mzBound, int smoothItrs, double smoothRTwindow, bool wavelet);
  getPepDataStrct(data, xmlFile, pngFileBase, cgiIndx, False, False, _ASAPRATIO_MZBOUND_, 10, 0.5, false); 

} 

void getPepDataStrct(pepDataStrct *data, char *xmlFile, char *pngFileBase, int cgiIndx, 
		     Boolean quantHighBG, Boolean zeroBG, double mzBound, bool wavelet) {
  void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		       char *pngFileBase, int cgiIndx, 
		       Boolean quantHighBG, Boolean zeroBG, 
		       double mzBound, int smoothItrs, double smoothRTwindow, bool wavelet);
  getPepDataStrct(data, xmlFile, pngFileBase, cgiIndx, quantHighBG, zeroBG, mzBound, 10, 0.5, wavelet); 

}


void getPepDataStrct(pepDataStrct *data, char *xmlFile, char *pngFileBase, int cgiIndx, 
		     Boolean quantHighBG, Boolean zeroBG, double mzBound, int smoothItrs, double smoothRTwindow, bool wavelet)
{
  lcSpectStrct lcSpect;
  int lcSize, cidScan;
  int cidIndx = data->cidIndx;
  int count, qState;
  double value[_ASAPRATIO_MXQ_];
  double error[_ASAPRATIO_MXQ_];
  double weight[_ASAPRATIO_MXQ_];
  int outliers[_ASAPRATIO_MXQ_];
  double timeWd = 3.; // time range for LC spectra in minutes
  int scanRange = 50;
  int pkIndx[2] = {0, 0};
  int valleyScan;
  double area;

  int i, j;

  if(data->indx < 0)
    return;

  if (mzBound <= 0) {
    mzBound =  _ASAPRATIO_MZBOUND_;
  }

  // get LC spectrum
  if(getPeakSpectra(&lcSpect, data, timeWd, xmlFile, mzBound, smoothItrs, smoothRTwindow, wavelet) == 0){
    data->indx = -1;
    return;
  }

  lcSize = lcSpect.rawSpectrum[data->chrg-1][cidIndx].size;
  cidScan = 0;
  while(cidScan < lcSize-1 
	&& lcSpect.expScanNums[cidScan] < data->scan) 
    ++cidScan;

  // convert valleys into scan numbers in LC/MS spectrum
  if(data->indx > 0) {
    for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
      for (j = 0; j < 2; ++j) {
	if(data->peaks[i][j].indx >= 0) {
	  if(lcSpect.expScanNums[0] > data->peaks[i][j].valley[1]
	     || lcSpect.expScanNums[lcSize-1] < data->peaks[i][j].valley[0]) {
	    data->peaks[i][j].peak = lcSize/2;
	    data->peaks[i][j].valley[0] = lcSize/2;
	    data->peaks[i][j].valley[1] = lcSize/2;
	  }
	  else {
	    //DDS: if the peak doesn't fall between the valleys define a new peak
	    if (data->peaks[i][j].peak < data->peaks[i][j].valley[0] ||
		data->peaks[i][j].peak > data->peaks[i][j].valley[1]) {
	      data->peaks[i][j].peak = (data->peaks[i][j].valley[0] + data->peaks[i][j].valley[1]) / 2;
	    }
	    valleyScan = 0;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].valley[0]) 
	      ++valleyScan;
	    data->peaks[i][j].valley[0] = valleyScan;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].peak) 
	      ++valleyScan;
	    data->peaks[i][j].peak = valleyScan;
	    while(valleyScan < lcSize-1
		  && lcSpect.expScanNums[valleyScan] 
		  < data->peaks[i][j].valley[1]) 
	      ++valleyScan;
	    data->peaks[i][j].valley[1] = valleyScan;
	    if(data->peaks[i][j].valley[0] >= data->peaks[i][j].valley[1]) {
	      data->peaks[i][j].peak = lcSize/2;
	      data->peaks[i][j].valley[0] = lcSize/2;
	      data->peaks[i][j].valley[1] = lcSize/2;
	    }
	  } // else {
	} //if(data->peaks[i][j].indx >= 0) {
      } //for (j = 0; j < 2; ++j) {
    } //for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
  } // if(data->indx > 0) {


  // get peak scan
  if(data->indx == 0) {
    // identified peptide 
    pkIndx[cidIndx] = -1;
    qState = data->chrg;
    if(data->peaks[qState-1][cidIndx].indx >= 0) {
      getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
		     lcSpect.rawSpectrum[qState-1][cidIndx],
		     lcSpect.fitSpectrum[qState-1][cidIndx],
		     cidScan, scanRange, data->areaFlag, quantHighBG, zeroBG);
      pkIndx[cidIndx] = qState - 1;
    } // if(data->peaks[qState-1][cidIndx].indx >= 0) {
    else {
      for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
	if(data->chrg > _ASAPRATIO_MXQ_/2) {
	  qState = data->chrg - i;
	  if(qState <= 0) qState += _ASAPRATIO_MXQ_;
	}	    
	else {
	  qState = data->chrg + i;
	  if(qState > _ASAPRATIO_MXQ_) qState -= _ASAPRATIO_MXQ_;
	}
	if(data->peaks[qState-1][cidIndx].indx >= 0) {
	  getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
			 lcSpect.rawSpectrum[qState-1][cidIndx],
			 lcSpect.fitSpectrum[qState-1][cidIndx],
			 cidScan, scanRange, data->areaFlag, quantHighBG, zeroBG);
	  if(data->peaks[qState-1][cidIndx].indx > 1) {
	    pkIndx[cidIndx] = qState - 1;
	    break;
	  }
	} // if(data->peaks[qState-1][cidIndx].indx >= 0) {
      } // for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
      if(pkIndx[cidIndx] < 0) {
	area = 0.;
	for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
	  if(data->peaks[i][cidIndx].indx > 0
	     && data->peaks[i][cidIndx].area[0] > area) {
	    area = data->peaks[i][cidIndx].area[0];
	    pkIndx[cidIndx] = i;
	  }
	}
      }
      if(pkIndx[cidIndx] < 0) 
	pkIndx[cidIndx] = data->chrg - 1;

      for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
	if(data->peaks[i][cidIndx].indx > 0
	   && pkIndx[cidIndx] != i) {
	  data->peaks[i][cidIndx].indx = 0;
	}
      }
    } // else {

    // peptide partner
    if(cidIndx == 0) { // light
      cidScan = data->peaks[pkIndx[cidIndx]][cidIndx].peak 
	+ data->eltn*(data->peaks[pkIndx[cidIndx]][cidIndx].valley[1]
		      -data->peaks[pkIndx[cidIndx]][cidIndx].valley[0])/4;
    }
    else { // heavy
      cidScan = data->peaks[pkIndx[cidIndx]][cidIndx].peak 
	- data->eltn*(data->peaks[pkIndx[cidIndx]][cidIndx].valley[1]
		      -data->peaks[pkIndx[cidIndx]][cidIndx].valley[0])/4;
    }
    cidScan = cidScan > 0 ? cidScan : 0;
    cidScan = cidScan < lcSize-1 ? cidScan : lcSize-1;

    cidIndx = (cidIndx+1)%2;
    pkIndx[cidIndx] = -1;
    qState = pkIndx[(cidIndx+1)%2] + 1;
    if(data->peaks[qState-1][cidIndx].indx >= 0) {
      getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
		     lcSpect.rawSpectrum[qState-1][cidIndx],
		     lcSpect.fitSpectrum[qState-1][cidIndx],
		     cidScan, scanRange, data->areaFlag, quantHighBG, zeroBG);
      if(data->peaks[qState-1][cidIndx].indx > 1) {
	pkIndx[cidIndx] = qState - 1;
      }
    } // if(data->peaks[qState-1][cidIndx].indx >= 0) {
    if(pkIndx[cidIndx] < 0) {
      for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
	if(data->chrg > _ASAPRATIO_MXQ_/2) {
	  qState = data->chrg - i;
	  if(qState <= 0) qState += _ASAPRATIO_MXQ_;
	}	    
	else {
	  qState = data->chrg + i;
	  if(qState > _ASAPRATIO_MXQ_) qState -= _ASAPRATIO_MXQ_;
	}
	if(data->peaks[qState-1][cidIndx].indx >= 0) {
	  getLCPeakStrct(&(data->peaks[qState-1][cidIndx]), 
			 lcSpect.rawSpectrum[qState-1][cidIndx],
			 lcSpect.fitSpectrum[qState-1][cidIndx],
			 cidScan, scanRange, data->areaFlag,quantHighBG, zeroBG);
	  if(data->peaks[qState-1][cidIndx].indx > 1) {
	    pkIndx[cidIndx] = qState - 1;
	    break;
	  }
	} // if(data->peaks[qState-1][cidIndx].indx >= 0) {
      } // for (i = 1; i < _ASAPRATIO_MXQ_; ++i){
    } // if(pkIndx[cidIndx] < 0) {
    if(pkIndx[cidIndx] < 0) {
      area = 0.;
      for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
	if(data->peaks[i][cidIndx].indx > 0
	   && data->peaks[i][cidIndx].area[0] > area) {
	  area = data->peaks[i][cidIndx].area[0];
	  pkIndx[cidIndx] = i;
	}
      }
    }
    if(pkIndx[cidIndx] < 0) 
      pkIndx[cidIndx] = pkIndx[(cidIndx+1)%2];

    for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
      if(data->peaks[i][cidIndx].indx > 0
	 && pkIndx[cidIndx] != i) {
	data->peaks[i][cidIndx].indx = 0;
      }
    }
  } //   if(data->indx == 0) {


  // get individual ratios
  area = 0.;
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    if(data->indx == 0) 
      data->pkCount[i] = 1;
    for (j = 0; j < 2; ++j) {
      if(data->peaks[i][j].indx == 0) {
	getLCPeakStrct(&(data->peaks[i][j]), lcSpect.rawSpectrum[i][j], 
		       lcSpect.fitSpectrum[i][j], 
		       data->peaks[pkIndx[j]][j].peak, scanRange, data->areaFlag, quantHighBG, zeroBG);
      }// if(data->peaks[i][j].indx == 0) 
      else if(data->peaks[i][j].indx > 0
	      && data->pkCount[i] == 1) { 
	data->peaks[i][j].indx 
	  = getLCSpectAreaAndTime(data->peaks[i][j].area, data->peaks[i][j].time,
				  lcSpect.rawSpectrum[i][j], lcSpect.fitSpectrum[i][j], 
				  data->peaks[i][j].peak, data->peaks[i][j].valley,
				  data->peaks[i][j].bckgrnd, data->areaFlag, quantHighBG);
      }
      if(data->indx == 0
	 && data->peaks[i][j].indx >= 0
	 && abs(data->peaks[i][j].peak-data->peaks[pkIndx[j]][j].peak)
	 > (data->peaks[i][j].valley[1]-data->peaks[i][j].valley[0])/4
	 + (data->peaks[pkIndx[j]][j].valley[1]-data->peaks[pkIndx[j]][j].valley[0])/4) 
	data->pkCount[i] = 0;
    } //for (j = 0; j < 2; ++j) {
    if(data->peaks[i][0].indx < 0 // at least 1 invalid peak
       || data->peaks[i][1].indx < 0
       || (data->peaks[i][0].indx == 1 // no valid area 
	   && data->peaks[i][1].indx == 1)) { 
      data->pkRatio[i] = -2;
      data->pkError[i] = 0.;      
      data->pkCount[i] = 0;      
    }
    else if(data->peaks[i][0].indx == 1) { // light invalid
      data->pkRatio[i] = 0.;
      data->pkError[i] = 0.;      
    }
    else if(data->peaks[i][1].indx == 1) { // heavy invalid
      data->pkRatio[i] = -1.;
      data->pkError[i] = 0.;
    }
    else { // both valid
      data->pkRatio[i] = data->peaks[i][0].area[0]/data->peaks[i][1].area[0];
      data->pkError[i] = data->pkRatio[i]
	*sqrt(data->peaks[i][0].area[1]*data->peaks[i][0].area[1]
	      /data->peaks[i][0].area[0]/data->peaks[i][0].area[0]
	      + data->peaks[i][1].area[1]*data->peaks[i][1].area[1]
	      /data->peaks[i][1].area[0]/data->peaks[i][1].area[0]);
    }
    if(data->pkCount[i] == 1
       && area < data->peaks[i][0].area[0]+data->peaks[i][1].area[0])
      area = data->peaks[i][0].area[0] + data->peaks[i][1].area[0];      
  } //for (i = 0; i < _ASAPRATIO_MXQ_; ++i){

  // strong data only
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    if (data->indx == 0
	&& data->pkCount[i] == 1
	&& data->peaks[i][0].area[0]+data->peaks[i][1].area[0] < 0.1*area)
      data->pkCount[i] = 0;
  }


  // get peptide output

  // collect valid date
  count = 0;
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    if (data->pkCount[i] == 1) {  
      value[count] = data->pkRatio[i];
      error[count] = data->pkError[i];
      weight[count] = data->peaks[i][0].area[0] 
	+ data->peaks[i][1].area[0];
      outliers[count] = 0; // assume not an outlier
      ++count;
    }
  }

  // calculate ratio and error
  if(count > 0) {
    // get ratio and error
    if(data->indx == 0) 
      getDataRatio(&(data->pepRatio[0]), &(data->pepRatio[1]), 
		   &(data->pepH2LRatio[0]), &(data->pepH2LRatio[1]), 
		   _ASAPRATIO_CONFL_,
		   value, error, weight, outliers, count, 1);
    else
      getDataRatio(&(data->pepRatio[0]), &(data->pepRatio[1]),  
		   &(data->pepH2LRatio[0]), &(data->pepH2LRatio[1]), 
		   _ASAPRATIO_CONFL_,
		   value, error, weight, outliers, count, 0);

    // set pepCount
    count = 0;
    for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
      if (data->pkCount[i] == 1) {  
	data->pkCount[i] = 1 - outliers[count];            
	++count;
      }
    }

    //get time and error
    for (i = 0; i < 2; ++i) {
      area = 0.;
      qState = data->chrg-1;      
      for (j = 0; j < _ASAPRATIO_MXQ_; ++j){
	if(data->peaks[j][i].indx == 2
	   && data->pkCount[j] == 1
	   && data->peaks[j][i].area[0] > area) {
	  qState = j;
	  area = data->peaks[j][i].area[0];
	}
      }
      data->pepTime[i][0] = data->peaks[qState][i].time[0];
      data->pepTime[i][1] = data->peaks[qState][i].time[1];
    }// for (i = 0; i < 2; ++i) {

    // get area
    data->pepArea = 0.;
    for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
      if (data->pkCount[i] == 1 &&
	  data->pepArea < data->peaks[i][0].area[0] + data->peaks[i][1].area[0])
	data->pepArea = data->peaks[i][0].area[0] + data->peaks[i][1].area[0];
    }
  } //  if(count > 0) {
  else {
    data->pepRatio[0] = -2.;
    data->pepRatio[1] = 0.;
    data->pepTime[0][0] = -1.;
    data->pepTime[0][1] = -1.;
    data->pepTime[1][0] = -1.;
    data->pepTime[1][1] = -1.;
    data->pepArea = -1.;
  }


  // generate .pngFiles
  if(cgiIndx == 1) {
    for (i = 0; i < _ASAPRATIO_MXQ_; ++i) {
      if(data->peaks[i][0].indx > 0 ||
	 data->peaks[i][1].indx > 0 ) {
	plotPepPeak(pngFileBase, *data, lcSpect, i+1);
      }
    }
  }
  // (assume == 2) generate js variables for use in flot.js
  else {
    for (i = 0; i < _ASAPRATIO_MXQ_; ++i)
      if(data->peaks[i][0].indx > 0 ||
	 data->peaks[i][1].indx > 0 )
	jsPlotPepPeak(*data, lcSpect, i+1);
  }


  // convert valleys into scan numbers in .dat file
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    for (j = 0; j < 2; ++j) {
      if(data->peaks[i][j].indx >= 0) {
	data->peaks[i][j].peak      = lcSpect.expScanNums[data->peaks[i][j].peak];
	data->peaks[i][j].valley[0] = lcSpect.expScanNums[data->peaks[i][j].valley[0]];
	data->peaks[i][j].valley[1] = lcSpect.expScanNums[data->peaks[i][j].valley[1]];
      }
    }
  }

  // reset data->indx;
  if(data->indx == 0)
    data->indx = 1;

  return;
}
  

// This function collects LC spectra for a pepDataStrct.  On success, it returns 1; on failure, 0.
int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, double timeWd, char *xmlFile) {
  return getPeakSpectra(lcSpect, data, timeWd, xmlFile, _ASAPRATIO_MZBOUND_, 10, 0.5, false);
}

int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, double timeWd, char *xmlFile,
		   double mzBound, int smoothItrs, double smoothRTwindow) {
  return getPeakSpectra(lcSpect, data, timeWd, xmlFile, _ASAPRATIO_MZBOUND_, smoothItrs, smoothRTwindow, false);
}

int getPeakSpectra(lcSpectStrct *lcSpect, pepDataStrct *data, double timeWd, char *xmlFile,
		   double mzBound, int smoothItrs, double smoothRTwindow, bool wavelet)
{
  xmlIndxStrct *getXmlIndx(char *xmlFile);
  void freeXmlIndx(xmlIndxStrct *xmlIndx);
  int getScanNum(double time, xmlIndxStrct *xmlIndx);
  spectStrct *getCombLCSpect(int firstScan, int lastScan, 
			     double *mzArray, int arraySize, double dMz, 
			     xmlIndxStrct *xmlIndx);

  xmlIndxStrct *xmlIndx;
  struct ScanHeaderStruct header;
  spectStrct *tmpSpect;
  double scanTime;
  long startScan, endScan;
  double mzArray[_ASAPRATIO_ISONUM_];
  double dMz;
  double mass[2];
  int lcScan, qState, size;
  
  long scan;

  int i, j, k;

  // initial steps
  if((xmlIndx = getXmlIndx(xmlFile)) == NULL){
    printf("Error In Reading %s.\n", xmlFile);
    fflush(stdout);
    return 0;
  }

  // get scanTime
  readHeader(xmlIndx->file, xmlIndx->scanIndx[data->scan], &header);
  scanTime = header.retentionTime/60.;

  // find scan range
  startScan = getScanNum(scanTime-timeWd, xmlIndx);
  endScan = getScanNum(scanTime+timeWd, xmlIndx);

  // find mass range
  mass[0] = data->msLight;
  mass[1] = data->msHeavy;

  
  // get spectra
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    qState = i + 1;
    for (j = 0; j < 2; ++j) {
      if(data->indx != 0 && data->peaks[i][j].indx == -1)
	continue;
      
      // get mzArray
      for (k = 0; k < _ASAPRATIO_ISONUM_; ++k) { 
	mzArray[k] = (mass[j]+k*_ASAPRATIO_M_)/qState + _ASAPRATIO_HM_*(qState - 1)/qState;
      }
      //      dMz = _ASAPRATIO_MZBOUND_ < _ASAPRATIO_M_/qState/2. ? 
      //	_ASAPRATIO_MZBOUND_ : _ASAPRATIO_M_/qState/2.;

      dMz = mzBound < _ASAPRATIO_M_/qState/2. ? 
      	mzBound : _ASAPRATIO_M_/qState/2.;
      

      // get raw LC spectrum
      if((tmpSpect = getCombLCSpect(startScan, endScan, mzArray, 
				    _ASAPRATIO_ISONUM_, dMz, xmlIndx)) == NULL) {
	data->peaks[i][j].indx = -1;
	continue;
      }
      else if(data->indx == 0)
	data->peaks[i][j].indx = 0;

      // get fitting spectrum
      lcSpect->rawSpectrum[i][j] = *tmpSpect;
      lcSpect->fitSpectrum[i][j] = *tmpSpect;
      lcSpect->fitSpectrum[i][j].smoothSpectFlx( smoothRTwindow, smoothItrs, 10,wavelet);
      delete tmpSpect;
    }//for (j = 0; j < 2; ++j) {
  } //  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
  

  // get expScanNums in mzXML file and cidScan in spectrum
  size = endScan - startScan + 1;
  lcSpect->expScanNums = (long *) calloc(size, sizeof(long));
  lcScan = 0;
  for (scan = startScan; scan <= endScan; ++scan) {
    readHeader(xmlIndx->file, xmlIndx->scanIndx[scan], &header);
    if (header.msLevel == 1 && 
	strcasecmp(header.scanType, "Zoom") != 0) {
      lcSpect->expScanNums[lcScan] = scan;
      ++lcScan;
    }
  }

  freeXmlIndx(xmlIndx);

  return 1;
}


// This function opens an xmlFile, gets scanIndx, and total scan number.
xmlIndxStrct *getXmlIndx(char *xmlFile)
{
  xmlIndxStrct *xmlIndx;
  ramp_fileoffset_t indexOffset;

  xmlIndx = (xmlIndxStrct *) calloc(1, sizeof(xmlIndxStrct));

  // open file
  if ((xmlIndx->file = rampOpenFile(xmlFile)) == NULL) {
    printf ("Could not open input file %s\n", xmlFile);
    fflush(stdout);
    free(xmlIndx);
    return NULL;
  }
  
  // Read the offset of the index
  indexOffset = getIndexOffset (xmlIndx->file);
  
  // Read the scan index into a vector, get LastScan
  xmlIndx->scanIndx = readIndex (xmlIndx->file, indexOffset, &(xmlIndx->totScan));
  
  return xmlIndx;
}


// This function frees a xmlIndxStrct.
void freeXmlIndx(xmlIndxStrct *xmlIndx)
{
  free(xmlIndx->scanIndx);
  rampCloseFile(xmlIndx->file);

  return;
}

// This function gets LC spectrum by summing up ion intensities in MS spectra within a series of mass windows.
spectStrct *getCombLCSpect(int firstScan, int lastScan, 
			   double *mzArray, int arraySize, double dMz, 
			   xmlIndxStrct *xmlIndx)
{
  spectStrct *spectrum;
  struct ScanHeaderStruct scanHeader;
  int iSt, iNd;
  int size;
  int dataNum;
  RAMPREAL *pPeaks;
  int peakCount;

  int i, j;

  // check on scan
  iSt = firstScan > 1 ? firstScan : 1;
  iNd = lastScan < xmlIndx->totScan ? lastScan : xmlIndx->totScan;
  if(iSt > iNd)
    return NULL;

  // get LC spectrum
  spectrum = new spectStrct(_MZXML_READER_MXSCANNUM_);
  size = 0;
  dataNum = 0;
  for (i = iSt; i <= iNd; ++i) {

    readHeader (xmlIndx->file, xmlIndx->scanIndx[i], &scanHeader);
    if (scanHeader.msLevel == 1 && 
	strcasecmp(scanHeader.scanType, "Zoom") != 0) { // LC/MS scan

      spectrum->xval[size] = scanHeader.retentionTime/60.;  // in minute
      spectrum->yval[size] = 0.;

      peakCount = 0;
      pPeaks = readPeaks(xmlIndx->file, xmlIndx->scanIndx[i]);
      
      if(pPeaks == NULL) {
	return NULL;
      }

      int l = 0;
      int r = scanHeader.peaksCount;
      // find monoisotopic peak
      RAMPREAL fMass;
      while (l < r) {
	int m = (l + r) / 2;
	fMass = pPeaks[m*2];
	if (fMass == mzArray[0]-dMz) { //probably an exact match is rare
	  l = r = m;
	  break;
	}
	else if (fMass > mzArray[0]-dMz) {
	  if (r == m) {
	    r = m - 1;
	  }
	  else {
	    r = m ;
	  }
	}
	else {
	  if (l == m) {
	    l = m + 1;
	  }
	  else {
	    l = m ;
	  }
	}

      }
      int peakIx = r * 2; //We should have the smallest m/z in the spectrum matching to the target range

      
      if(peakIx <= scanHeader.peaksCount*2 && pPeaks[peakIx] >= mzArray[0]-dMz) {
	// Sum the signal over all peaks falling within the target range around each isotope
	j = 0;
	while (j < arraySize && pPeaks[peakIx] != -1 && pPeaks[peakIx] <= mzArray[arraySize-1]+dMz ) {
	  if (pPeaks[peakIx] <= mzArray[j]+dMz
	       && pPeaks[peakIx] >= mzArray[j]-dMz) {
	      spectrum->yval[size] += pPeaks[peakIx+1];
	    peakIx += 2;
	    }
	  else if (pPeaks[peakIx] > mzArray[j]+dMz) {
	    j++;
	  }
	  else if (pPeaks[peakIx] < mzArray[j]-dMz) {
	  peakIx += 2;
	}
	  else {
	    peakIx += 2; // Should never get here
	  }

	}
      }

      //free(pPeaks);
      
      if(size < spectrum->size && spectrum->yval[size] > 0.)
	++dataNum;
      ++size;
    } //     if (scanHeader.msLevel == 1) { // LC/MS scan

  } // for (i = iSt; i <= iNd; ++i) {

  if(dataNum < 1) {
    delete spectrum;
    return NULL;
  }
  else {
    spectrum->size = size;
    spectrum->xval = (double *) realloc(spectrum->xval, spectrum->size*sizeof(double));
    spectrum->yval = (double *) realloc(spectrum->yval, spectrum->size*sizeof(double));
    
    return spectrum;
  }
}

// For a given time, this function returns the closest scan number with fabs(scanTime-time) minimium.
int getScanNum(double time, xmlIndxStrct *xmlIndx)
{
  struct ScanHeaderStruct header;
  int scanNum;
  int bnd[2];
  double scanTime[2];
  double tmpTime;

  // 1st spect
  bnd[0] = 1;
  readHeader(xmlIndx->file, xmlIndx->scanIndx[bnd[0]], &header);
  scanTime[0] = header.retentionTime/60.;
  if(scanTime[0] >= time)
    return bnd[0];

  // last spect
  bnd[1] = xmlIndx->totScan;
  readHeader(xmlIndx->file, xmlIndx->scanIndx[bnd[1]], &header);
  scanTime[1] = header.retentionTime/60.;
  if(scanTime[1] <= time)
    return bnd[1];

  // search for scanNum
  while(bnd[1]-bnd[0] > 1) {
    scanNum = (bnd[0]+bnd[1])/2;
    readHeader(xmlIndx->file, xmlIndx->scanIndx[scanNum], &header);
    tmpTime = header.retentionTime/60.;
    if(tmpTime == time)
      return scanNum;
    else if(tmpTime < time) {
      bnd[0] = scanNum;
      scanTime[0] = tmpTime;
    }
    else {
      bnd[1] = scanNum;
      scanTime[1] = tmpTime;
    }
  }

  if(fabs(scanTime[0]-time) <= fabs(scanTime[1]-time))
    return bnd[0];
  else
    return bnd[1];
}


// This function gets a lcPeakStrct for a spectrum
void getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, int scan, int scanRange) {
  getLCPeakStrct(peak,  rawSpectrum, fitSpectrum,  scan,  scanRange, 0, False, False);
}

void getLCPeakStrct(lcPeakStrct *peak, const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, int scan, int scanRange,
		    int areaType, Boolean quantHighBG, Boolean zeroBG)
{
  int ranges[2];

  // get peakScan and valley
  getLCSpectPeakAndValleys(&(peak->peak), peak->valley, fitSpectrum, scan, 0.); 
  
  // get background
  ranges[0] = peak->valley[0] - scanRange;
  ranges[0] = ranges[0] > 0 ? ranges[0] : 0;
  ranges[1] = peak->valley[1] + scanRange;
  ranges[1] = ranges[1] < rawSpectrum.size-1 ? ranges[1] : rawSpectrum.size-1;
  peak->bckgrnd = getLCSpectBackground(rawSpectrum, fitSpectrum, ranges, peak->valley, zeroBG);
  
  // reset valleys
  getLCSpectPeakAndValleys(&(peak->peak), peak->valley, fitSpectrum, 
			   peak->peak, peak->bckgrnd); 
  
  // get area and time
  peak->indx = getLCSpectAreaAndTime(peak->area, peak->time, 
				     rawSpectrum, fitSpectrum, 
				     peak->peak, peak->valley, 
				     peak->bckgrnd, areaType, quantHighBG);

  return;
}


// This function finds the peak and valleys of a LC spectrum.
void getLCSpectPeakAndValleys(int *peakScan, int *valleyScan, const spectStrct &spectrum, int scan, double background) 
{
  int i;
  
  // find peak scan
  *peakScan = spectPeak(spectrum, scan, -1);

  // find valley positions
  valleyScan[0] = findNextValley(spectrum, *peakScan, -1);
  if (valleyScan[0] == -1) { // if no valley found, set to boundary
    valleyScan[0] = 0;
  }
  valleyScan[1] = findNextValley(spectrum, *peakScan, +1);
  if (valleyScan[1] == -1) { // if no valley found, set to boundary
    valleyScan[1] = spectrum.size-1;
  }

  // find positions at background level
  for (i = *peakScan; i >= valleyScan[0]; --i){
    if(spectrum.yval[i] < background) {
      break;
    }
  }
  valleyScan[0] = i + 1;
  for (i = *peakScan; i <= valleyScan[1]; ++i){
    if(spectrum.yval[i] < background)
      break;
  }
  valleyScan[1] = i - 1;

  // check valley
  if(valleyScan[0] > valleyScan[1]) {
    valleyScan[0] = *peakScan;
    valleyScan[1] = *peakScan;
  }

  return;
}


// This function finds the peak at "position" in a spectrum.  
// If "position" is a valley, then find the left peak (when "direction = -1") or the right one (when "direction = 1").
int spectPeak(const spectStrct &spectrum, int position, int direction)
{
  int right, left;
  int i;

  // determine slope on the right hand side 
  i = 1;
  while(position+i < spectrum.size &&
	spectrum.yval[position] == spectrum.yval[position+i]){
    ++i;
  } 
  if (position+i >= spectrum.size) {
    right = 0;
  }
  else if (spectrum.yval[position] < spectrum.yval[position+i]){
    right = 1;
  }
  else
    right = -1;

  // determine slope on the left hand side 
  i = -1;
  while(position+i >= 0 &&
	spectrum.yval[position] == spectrum.yval[position+i]){
    --i;
  } 
  if (position+i < 0) {
    left = 0;
  }
  else if (spectrum.yval[position] < spectrum.yval[position+i]){
    left = 1;
  }
  else
    left = -1;

  // determine "direction"
  if((right == -1 && left == -1) || (right == 0 && left == 0)) {
    return position;
  }
  else if (right == -1 && left == 0) {
    return 0;
  }
  else if (right == 0 && left == -1) {
    return spectrum.size-1;
  }
  else if (right <= 0 && left == 1) {
    direction = -1;
  }
  else if (right == 1 && left <= 0) {
    direction = 1;
  }

  // search for peak
  while (position+direction >= 0 && position+direction < spectrum.size
	 && spectrum.yval[position] <= spectrum.yval[position+direction]) {
    position += direction;
  }

  // return index for peak
  if (position+direction < 0) 
    return 0;
  else if (position+direction >= spectrum.size) 
    return spectrum.size-1;
  else 
    return position;
}

// This function finds the next valley from left (when "direction = -1") a
// or right (when "direction = 1"), starting at "position" in a spectrum.  
int findNextValley(const spectStrct &spectrum, int position, int direction)
{
  position += direction;

  // search for valley
  while (position-direction >= 0 && position-direction < spectrum.size 
	 && position+direction >= 0 && position+direction < spectrum.size 
	 && (spectrum.yval[position] > spectrum.yval[position-direction] ||
	     spectrum.yval[position] >= spectrum.yval[position+direction])) {
    position += direction;
  }

  // return index for valley
  if (position-direction < 0 || position-direction >= spectrum.size 
      || position+direction < 0 || position+direction >= spectrum.size) {
    return -1;
  } 
  else if(spectrum.yval[position] <= spectrum.yval[position-direction] 
	  && spectrum.yval[position] < spectrum.yval[position+direction]){
    return position;
  }
  else {
    return -1;
  }
}


// This function gets the background of a LC spectrum.
double getLCSpectBackground(const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, int *ranges, int *valleyScan) {
 return  getLCSpectBackground( rawSpectrum, fitSpectrum, ranges, valleyScan, False);
}

double getLCSpectBackground(const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, int *ranges, int *valleyScan,
			    Boolean zeroBG)
{
  double noise, background, ave;
  int count, count2;
  int i;

  // check on ranges
  ranges[0] = ranges[0] > 0 ? ranges[0] : 0;
  ranges[1] = ranges[1] < rawSpectrum.size-1 ?
    ranges[1] : rawSpectrum.size-1;

  // get noise level
  noise = 0.;
  count = 0;
  for (i = ranges[0]; i <= ranges[1]; ++i){
    if (rawSpectrum.yval[i] > 0.) {
      noise += (rawSpectrum.yval[i]-fitSpectrum.yval[i])
	*(rawSpectrum.yval[i]-fitSpectrum.yval[i]);
      ++count;
    }
  }
  if(count > 0)
    noise = sqrt(noise/count);
  else
    noise = 1.;

  // get background
  background = 0.;
  ave = noise;
  while(fabs(background-ave) > 0.01*noise) {
    background = ave;
    ave = 0.;
    count = 0;
    count2 = 0;
    for (i = ranges[0]; i <= ranges[1]; ++i){
      if(i < valleyScan[0] || i > valleyScan[1]) {
	if (rawSpectrum.yval[i] > 0.
	    && rawSpectrum.yval[i] < background + noise) {
	  ave += rawSpectrum.yval[i];
	  ++count;
	}
	else if (rawSpectrum.yval[i] >= background + noise) 
	  ++count2;
      }
    }// for (i = ranges[0]; i <= ranges[1]; ++i){
    if(count > 0)
      ave /= count;
    else if(count2 > 0)
      ave = background + 0.2*noise;
    else
      ave = 0.;
  }

  if (zeroBG) {
    background = 0.;
  }
  else {
  background = ave;
  }

  return background;
}


// This function finds the area and time of a LC peak. It returns 2 if the peak passes certain tests, 1 if not.
int getLCSpectAreaAndTime(double *area, double *peakTime, const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, 
			  int pkScan, int *valleyScan, double background) {
  return getLCSpectAreaAndTime(area, peakTime, rawSpectrum, fitSpectrum, pkScan, valleyScan,  background, 0, False);

 }
int getLCSpectAreaAndTime(double *area, double *peakTime, const spectStrct &rawSpectrum, const spectStrct &fitSpectrum, 
			  int pkScan, int *valleyScan, double background, int areaType, Boolean quantHighBG)
{
  int areaIndx = 2;
  double rArea, fArea;
  double rVal, fVal;
  double mxVal, ave, areaErr;
  int pkTime1, pkTime2;
  int count;
  int i;
  
  // get mxVal
  mxVal = 0.;
  for (i = valleyScan[0]; i <= valleyScan[1]; ++i) {
    mxVal = mxVal > fitSpectrum.yval[i] ?
      mxVal : fitSpectrum.yval[i];
  }

  // get peak area and error
  rArea = 0.;
  fArea = 0.;
  areaErr = 0.;
  for (i = valleyScan[0]; i <= valleyScan[1]; ++i) {
    rVal = 0.5*(rawSpectrum.yval[i+1]+rawSpectrum.yval[i]) - background;
    rVal = rVal > 0. ? rVal : 0.;
    rArea += rVal*(rawSpectrum.xval[i+1]-rawSpectrum.xval[i]);
    fVal = 0.5*(fitSpectrum.yval[i+1]+fitSpectrum.yval[i]) - background;
    fVal = fVal > 0. ? fVal : 0.;
    fArea += fVal*(fitSpectrum.xval[i+1]-fitSpectrum.xval[i]);
    areaErr += 0.5*(rVal-fVal)*(rVal-fVal)
      *(fitSpectrum.xval[i+1]-fitSpectrum.xval[i])
      *(fitSpectrum.xval[i+1]-fitSpectrum.xval[i]);
  }

  // check error
  if (areaType == 1) {
    area[0] = rArea;
  }
  else if (areaType == 2) {
    area[0] = fArea;
  }
  else {
  area[0] = (fArea+rArea)/2.;
  }

  area[1] = sqrt(areaErr);
  if(rArea <= 0. || fArea <= 0.) {
    areaIndx = 1;
  }
  else if (quantHighBG) {
    areaIndx = 2;
  }
  else if(area[0] < area[1] || mxVal < 2.*background) {
    areaIndx = 1;
  }

  // find peakTime 
  peakTime[0] = fitSpectrum.xval[pkScan];
  
  // find error in peakTime
  pkTime1 = pkScan;
  while(pkTime1 >= valleyScan[0] // width at half peak height: left
	&& fitSpectrum.yval[pkTime1]-background 
	> 0.5*(fitSpectrum.yval[pkScan]-background)) {
    --pkTime1;
  }
  pkTime2 = pkScan;
  while(pkTime2 < valleyScan[1] // width at half peak height: right
	&& fitSpectrum.yval[pkTime2]-background 
	> 0.5*(fitSpectrum.yval[pkScan]-background)) {
    ++pkTime2;
  }
  if(pkTime2 <= pkTime1)
    peakTime[1] = 0.;
  else {
    peakTime[1] = 0.5*(fitSpectrum.xval[pkTime2]-fitSpectrum.xval[pkTime1]);
  }

  // double check on bad data
  if(area[0] > 0. && mxVal < 5.*background) {
    count = 0;
    for (i = pkTime1; i <= pkTime2; ++i){
      if(rawSpectrum.yval[i] < background) {
	++count;
      }
    }
    if(count > (pkTime2-pkTime1+1)/4) {
      count = 0;
      ave = fitSpectrum.yval[pkTime1] < fitSpectrum.yval[pkTime2]?
	fitSpectrum.yval[pkTime1] : fitSpectrum.yval[pkTime2];
      for (i = pkTime1; i <= pkTime2; ++i)
	if(rawSpectrum.yval[i] > ave)
	  ++count;
      if (quantHighBG) {
	areaIndx = 2;
      }
      else if(count < (pkTime2-pkTime1+1)/2) {
	areaIndx = 1;
      }
    }
  } // if(area[0] > 0. && mxVal < 5.*background) {

  return areaIndx;
}


// This function generates pngFiles for LC/MS spectra.
void plotPepPeak(const char *pngFileBaseIn, const pepDataStrct &data, 
		 const lcSpectStrct &lcSpect, int qState)
{
  FILE *file, *file2;

  char *pngFile;
  char *gnuFile;
  char *datFiles[3][2];
  //char *gnuCommand;
  int plotRange[2];
  long expRange[2];
  double xtics = 0.5;
  int cidIndx = data.cidIndx;
  int size = lcSpect.rawSpectrum[data.chrg-1][cidIndx].size;
  long cidTime = data.scan;
  double mxPeak, yMx, yOrder;
  int lcScan;

  int i, j;

  char *pngFileBase = strdup(pngFileBaseIn);
  fixPath(pngFileBase,0); // pretty up the path separators etc

  // get timeRange
  plotRange[0] = data.peaks[data.chrg-1][cidIndx].valley[0];
  plotRange[1] = data.peaks[data.chrg-1][cidIndx].valley[1];
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    for (j = 0; j < 2; ++j) {
      if(data.peaks[i][j].indx >= 0) {
	if(plotRange[0] > data.peaks[i][j].valley[0])
	  plotRange[0] = data.peaks[i][j].valley[0];
	if(plotRange[1] < data.peaks[i][j].valley[1])
	  plotRange[1] = data.peaks[i][j].valley[1];
      }
    }
  }
  lcScan = plotRange[0];
  while(lcScan > 0 
	&& lcSpect.expScanNums[lcScan] 
	> lcSpect.expScanNums[plotRange[0]]-_ASAPRATIO_EXPLCRANGE_/5)
    --lcScan;
  plotRange[0] = lcScan;
  lcScan = plotRange[1];
  while(lcScan < size-1
	&& lcSpect.expScanNums[lcScan] 
	< lcSpect.expScanNums[plotRange[1]]+_ASAPRATIO_EXPLCRANGE_/5)
    ++lcScan;
  plotRange[1] = lcScan;

  // get expRange
  expRange[0] = 5*(lcSpect.expScanNums[plotRange[0]]/5);
  expRange[1] = 5*(lcSpect.expScanNums[plotRange[1]]/5+1);

  // get datFiles
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 2; ++j) {
      datFiles[i][j] = (char *) calloc(strlen(pngFileBase)+20, sizeof(char));
      sprintf(datFiles[i][j], "%s_%d_%d_%d.dat", pngFileBase, qState, i, j);
    }
  }

  // get pngFile
  pngFile = (char *) calloc(strlen(pngFileBase)+10, sizeof(char));
  sprintf(pngFile, "%s_%d.png", pngFileBase, qState);

  // get gnuFile
  gnuFile = (char *) calloc(strlen(pngFileBase)+10, sizeof(char));
  sprintf(gnuFile, "%s_%d.gnu", pngFileBase, qState);
  if((file = fopen(gnuFile, "w")) == NULL) {
    printf("Cannot Write to File \"%s\"!\n", gnuFile);
    fflush(stdout);
    exit(0);
  }

  fprintf(file, "set term png small\n");
  fprintf(file, "set output \'%s\'\n", pngFile);
  fprintf(file, "set format y \"%%7.2g\"\n");
  fprintf(file, "set key samplen -1\n");
  fprintf(file, "set multiplot\n");
  for (j = 0; j < 2; ++j) { 
    if(data.peaks[qState-1][j].indx >= 0) {
      fprintf(file, "set size 1.,0.5\n");
      if(j == 0) { // light
	fprintf(file, "set origin 0,0.5\n");
      }
      else { // heavy
	fprintf(file, "set origin 0,0\n");
      }
      fprintf(file, "set xrange [%ld:%ld]\n", expRange[0], expRange[1]);
      fprintf(file, "set yrange [0:*]\n");

      // datFile: rawSpectrum
      yMx = 0.;
      mxPeak = 0.;
      if((file2 = fopen(datFiles[0][j], "w")) == NULL) {
	printf("Cannot Write to File \"%s\"!\n", datFiles[0][j]);
	fflush(stdout);
	exit(0);
      }
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
	   i < plotRange[1]+10 && i < lcSpect.rawSpectrum[qState-1][j].size; 
	   ++i) {
	if(i > data.peaks[qState-1][j].valley[0]
	   && i < data.peaks[qState-1][j].valley[1]
	   && mxPeak < lcSpect.rawSpectrum[qState-1][j].yval[i])
	  mxPeak = lcSpect.rawSpectrum[qState-1][j].yval[i];
	if(i > plotRange[0] && i < plotRange[1] 
	   && yMx < lcSpect.rawSpectrum[qState-1][j].yval[i])
	  yMx = lcSpect.rawSpectrum[qState-1][j].yval[i];
	fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], lcSpect.rawSpectrum[qState-1][j].yval[i]);
      }
      fclose(file2);

      if(yMx > 4.*mxPeak && mxPeak > 0.) {
	yOrder = pow(10., (int)log10(mxPeak));
	mxPeak = 2.*ceil(mxPeak/yOrder)*yOrder;
	fprintf(file, "set yrange [0:%f]\n", mxPeak);
      }

      // datFile: fitSpectrum
      if((file2 = fopen(datFiles[1][j], "w")) == NULL) {
	printf("Cannot Write to File \"%s\"!\n", datFiles[1][j]);	
	fflush(stdout);
	exit(0);
      }
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
	   i < plotRange[1]+10 && i < lcSpect.fitSpectrum[qState-1][j].size; 
	   ++i)
	fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], lcSpect.fitSpectrum[qState-1][j].yval[i]);
      fclose(file2);

      // datFile: peakSpectrum
      if(data.peaks[qState-1][j].indx > 1 
	 && data.pkCount[qState-1] == 1) {
	if((file2 = fopen(datFiles[2][j], "w")) == NULL) {
	  printf("Cannot Write to File \"%s\"!\n", datFiles[2][j]);
	  fflush(stdout);
	  exit(0);
	}
	for (i = data.peaks[qState-1][j].valley[0];
	     i <= data.peaks[qState-1][j].valley[1]; ++i) {
	  fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], data.peaks[qState-1][j].bckgrnd);

	  float inten;
	  if(data.areaFlag == 1)
	    if (lcSpect.rawSpectrum[qState-1][j].yval[i] > data.peaks[qState-1][j].bckgrnd)
	      inten = lcSpect.rawSpectrum[qState-1][j].yval[i];
	    else
	      inten = data.peaks[qState-1][j].bckgrnd;
	  else if(data.areaFlag == 2)
	    if (lcSpect.fitSpectrum[qState-1][j].yval[i] > data.peaks[qState-1][j].bckgrnd)
	      inten = lcSpect.fitSpectrum[qState-1][j].yval[i];
	    else
	      inten = data.peaks[qState-1][j].bckgrnd;
	  else
	    if (0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i] + lcSpect.fitSpectrum[qState-1][j].yval[i])
		> data.peaks[qState-1][j].bckgrnd)
	      inten = 0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i] + lcSpect.fitSpectrum[qState-1][j].yval[i]);
	    else
	      inten = data.peaks[qState-1][j].bckgrnd;

	  fprintf(file2, "%ld %f\n", lcSpect.expScanNums[i], inten);

	}
	fclose(file2);
      }

      // cid indx
      if(qState == data.chrg && j == cidIndx) {
	fprintf(file, "set arrow from  %ld, graph 0 to %ld, graph 0.9 head size screen 0.01,90 lc rgb '#ff5f00' lw 2\n",
		cidTime, cidTime);
	  fprintf(file, "set label \"CID\" at %ld, graph 0.95 tc rgb '#ff5f00'\n",
		cidTime);
      }
      else {
	fprintf(file, "set noarrow\n");
        fprintf(file, "set arrow from  %ld, graph 0 to %ld, graph 0.9 nohead lc rgb '#5e6a71' lw 1\n",
                cidTime, cidTime);
      	fprintf(file, "set nolabel\n");
      }

      // peakSpectrum
      if(data.peaks[qState-1][j].indx > 1
	 && data.pkCount[qState-1] == 1) 
	fprintf(file, "pl '%s't \"\"w l lc rgb '#22eceb' lw 2, ", datFiles[2][j]);
      else 
      	fprintf(file, "pl ");

      // rawSpectrum 
      double tmp_mz;
      if(j == 0) {
	tmp_mz = (data.msLight+(qState-1)*_ASAPRATIO_HM_)/(double)qState;
	fprintf(file, "'%s't \"light +%d, m/z %.4f\"w l lc rgb '#ff0000' lw 2, ",
		datFiles[0][j], qState, tmp_mz);
      }
      else { 
	tmp_mz= (data.msHeavy+(qState-1)*_ASAPRATIO_HM_)/(double)qState;
	fprintf(file, "'%s't \"heavy +%d, m/z %.4f\"w l lc rgb '#ff0000' lw 2, ",
		datFiles[0][j], qState, tmp_mz);
      }
      // fitSpectrum
      fprintf(file, "'%s't \"\"w l lc rgb '#006dff' lw 2, ", datFiles[1][j]);

      // background
      fprintf(file, "%f t \"\"w l lc rgb '#82003b' lw 2\n", 
	      data.peaks[qState-1][j].bckgrnd);
    }
  }

  fprintf(file, "set nomultipl\n");
  fprintf(file, "unset output\n");
  fprintf(file, "quit\n");  
  fflush(file);
  fclose(file);

  // get image
  std::string gnuCommand(GNUPLOT_BINARY);
  gnuCommand += " ";
  gnuCommand += gnuFile;

  verified_unlink(pngFile);
  verified_system(gnuCommand.c_str()); // system() with verbose error check	

  //Sleep Min 30 secs Max 30 secs
  FILE* test;
  int count = 0;
  while (count < 30 && (test = fopen(pngFile,"r"))==NULL) {
    sleep(1);
  }

  if (count >= 30) {
    fprintf(stderr,"WARNING: Max Timeout reached waiting for gnuplot ... ");
  }
  else if ( test!=NULL ) {
    fclose(test);
  }
	
  // clear
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 2; ++j) {
      verified_unlink(datFiles[i][j]);
      free(datFiles[i][j]);
    }
  }

  verified_unlink(gnuFile);
  free(gnuFile);

  free(pngFile);
  free(pngFileBase);

  return;
}





// Write javascript variables with peak data for display in flot.js
void jsPlotPepPeak(const pepDataStrct &data, const lcSpectStrct &lcSpect, int qState)
{
  int i, j;


  // /////////////////// need this?
  int cidIndx = data.cidIndx; //?
  int size = lcSpect.rawSpectrum[data.chrg-1][cidIndx].size; //?
  int lcScan;
  int plotRange[2];
  long expRange[2];

  // get timeRange
  plotRange[0] = data.peaks[data.chrg-1][cidIndx].valley[0];
  plotRange[1] = data.peaks[data.chrg-1][cidIndx].valley[1];
  for (i = 0; i < _ASAPRATIO_MXQ_; ++i){
    for (j = 0; j < 2; ++j) {
      if(data.peaks[i][j].indx >= 0) {
        if(plotRange[0] > data.peaks[i][j].valley[0])
          plotRange[0] = data.peaks[i][j].valley[0];
        if(plotRange[1] < data.peaks[i][j].valley[1])
          plotRange[1] = data.peaks[i][j].valley[1];
      }
    }
  }
  lcScan = plotRange[0];
  while(lcScan > 0
        && lcSpect.expScanNums[lcScan]
        > lcSpect.expScanNums[plotRange[0]]-_ASAPRATIO_EXPLCRANGE_/5)
    --lcScan;
  plotRange[0] = lcScan;
  lcScan = plotRange[1];
  while(lcScan < size-1
        && lcSpect.expScanNums[lcScan]
        < lcSpect.expScanNums[plotRange[1]]+_ASAPRATIO_EXPLCRANGE_/5)
    ++lcScan;
  plotRange[1] = lcScan;

  // get expRange
  expRange[0] = 5*(lcSpect.expScanNums[plotRange[0]]/5);
  expRange[1] = 5*(lcSpect.expScanNums[plotRange[1]]/5+1);
  // /////////////////// need this?


  printf("<script type='text/javascript'>\n");
  printf(" var raw_data_%d = [];\n",qState);
  printf(" var fit_data_%d = [];\n",qState);
  printf(" var peaks_data_%d = [];\n",qState);

  for (j = 0; j < 2; ++j) {
    if(data.peaks[qState-1][j].indx >= 0) {
      // rawSpectrum
      printf(" raw_data_%d[%d] = [ \n", qState,j);
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
           i < plotRange[1]+10 && i < lcSpect.rawSpectrum[qState-1][j].size;
           ++i)
	printf(" [%ld,%f],\n", lcSpect.expScanNums[i], lcSpect.rawSpectrum[qState-1][j].yval[i]);
      printf(" [] ];\n");

      // fitSpectrum
      printf(" fit_data_%d[%d] = [ \n", qState,j);
      for (i = plotRange[0]-10 > 0 ? plotRange[0]-10 : 0;
           i < plotRange[1]+10 && i < lcSpect.fitSpectrum[qState-1][j].size;
           ++i)
	printf("[%ld,%f],\n", lcSpect.expScanNums[i], lcSpect.fitSpectrum[qState-1][j].yval[i]);
      printf(" [] ];\n");

      // peakSpectrum
      if(data.peaks[qState-1][j].indx > 1
         && data.pkCount[qState-1] == 1) {
	printf(" peaks_data_%d[%d] = [ \n", qState,j);
        for (i = data.peaks[qState-1][j].valley[0];
             i <= data.peaks[qState-1][j].valley[1]; ++i) {
          printf("[%ld,%f],\n", lcSpect.expScanNums[i], data.peaks[qState-1][j].bckgrnd);

          float inten;
          if(data.areaFlag == 1)
            if (lcSpect.rawSpectrum[qState-1][j].yval[i] > data.peaks[qState-1][j].bckgrnd)
              inten = lcSpect.rawSpectrum[qState-1][j].yval[i];
            else
              inten = data.peaks[qState-1][j].bckgrnd;
          else if(data.areaFlag == 2)
            if (lcSpect.fitSpectrum[qState-1][j].yval[i] > data.peaks[qState-1][j].bckgrnd)
              inten = lcSpect.fitSpectrum[qState-1][j].yval[i];
            else
              inten = data.peaks[qState-1][j].bckgrnd;
          else
            if (0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i] + lcSpect.fitSpectrum[qState-1][j].yval[i])
                > data.peaks[qState-1][j].bckgrnd)
              inten = 0.5*(lcSpect.rawSpectrum[qState-1][j].yval[i] + lcSpect.fitSpectrum[qState-1][j].yval[i]);
            else
              inten = data.peaks[qState-1][j].bckgrnd;

	  printf("[%ld,%f],\n", lcSpect.expScanNums[i], inten);
        }
	printf(" [] ];\n");
      }
    }
  }

  printf("</script>\n");
}
