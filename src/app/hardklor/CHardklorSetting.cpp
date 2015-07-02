#include "CHardklorSetting.h"

using namespace std;

CHardklorSetting::CHardklorSetting(){
  noBase=false;
	iAnalysis=false;
  maxCharge=5;
	minCharge=1;
  depth=3;
  peptide=10;
  smooth=0;
  corr=0.85;
  sn=1.0;
  scan.iLower=0;
  scan.iUpper=0;
  window.dLower=0;
  window.dUpper=0;
  strcpy(inFile,"");
  strcpy(outFile,"");
	strcpy(MercuryFile,"");
  strcpy(HardklorFile,"");
	algorithm=FastFewestPeptides;
	variant = new vector<CHardklorVariant>;

	msType=FTICR;
	res400=100000;
	winSize=4.0;

	chargeMode='Q';
	snWindow=250.0;
	skipZero=true;
	noSplit=false;
	sl=2;
	fileFormat=dunno;
	mzXMLFilter=MS1;
	distArea=false;
	strcpy(formula,"");

	centroid=false;
	staticSN=true;
  xml=false;
	reducedOutput=false;

  //ppMatch=1;
  //ppWin=1;
  //noiseMatch=1;
  //noiseWindow=3;
  ppm=10.0;
  //sna=0;
	boxcar=0;
	boxcarFilter=0;

  //rawAvg=false;
  //rawAvgWidth=1;
  //rawAvgCutoff=1000;

  strcpy(rawFilter,"");
}

CHardklorSetting::CHardklorSetting(const CHardklorSetting& c){
  int i;

	//Copy variant list
  variant = new vector<CHardklorVariant>;
  for(i=0;i<c.variant->size();i++) variant->push_back(c.variant->at(i));
	
  //Copy other data memebers
  noBase=c.noBase;
	iAnalysis=c.iAnalysis;
  maxCharge=c.maxCharge;
	minCharge=c.minCharge;
  depth=c.depth;
  peptide=c.peptide;
  smooth=c.smooth;
  corr=c.corr;
  sn=c.sn;
  scan.iLower=c.scan.iLower;
  scan.iUpper=c.scan.iUpper;
  window.dLower=c.window.dLower;
  window.dUpper=c.window.dUpper;
  strcpy(inFile,c.inFile);
  strcpy(outFile,c.outFile);
  strcpy(MercuryFile,c.MercuryFile);
  strcpy(HardklorFile,c.HardklorFile);
	algorithm=c.algorithm;
	msType=c.msType;
	res400=c.res400;
	winSize=c.winSize;

	chargeMode=c.chargeMode;
	snWindow=c.snWindow;
	skipZero=c.skipZero;
	noSplit=c.noSplit;
	sl=c.sl;
	fileFormat=c.fileFormat;
	mzXMLFilter=c.mzXMLFilter;
	distArea=c.distArea;
	strcpy(formula,c.formula);

	centroid = c.centroid;
	staticSN = c.staticSN;
  xml = c.xml;
	reducedOutput = c.reducedOutput;

  //ppMatch=c.ppMatch;
  //ppWin=c.ppWin;
  //noiseMatch=c.noiseMatch;
  //noiseWindow=c.noiseWindow;
  ppm=c.ppm;
  //sna=c.sna;
	boxcar=c.boxcar;
	boxcarFilter=c.boxcarFilter;

  //rawAvg=c.rawAvg;
  //rawAvgWidth=c.rawAvgWidth;
  //rawAvgCutoff=c.rawAvgCutoff;
  strcpy(rawFilter,c.rawFilter);
}
  
CHardklorSetting::~CHardklorSetting(){
	delete variant;
}

CHardklorSetting& CHardklorSetting::operator=(const CHardklorSetting& c){
  int i;
  if (this!=&c){
		delete variant;
    variant = new vector<CHardklorVariant>;
    for(i=0;i<c.variant->size();i++){
      variant->push_back(c.variant->at(i));
    }
    noBase=c.noBase;
		iAnalysis=c.iAnalysis;
    maxCharge=c.maxCharge;
		minCharge=c.minCharge;
    depth=c.depth;
		peptide=c.peptide;
    smooth=c.smooth;
    corr=c.corr;
    sn=c.sn;
    scan.iLower=c.scan.iLower;
    scan.iUpper=c.scan.iUpper;
    window.dLower=c.window.dLower;
    window.dUpper=c.window.dUpper;
    strcpy(inFile,c.inFile);
    strcpy(outFile,c.outFile);
    strcpy(MercuryFile,c.MercuryFile);
    strcpy(HardklorFile,c.HardklorFile);
		algorithm=c.algorithm;
		msType=c.msType;
		res400=c.res400;
		winSize=c.winSize;

		chargeMode=c.chargeMode;
		snWindow=c.snWindow;
		skipZero=c.skipZero;
		noSplit=c.noSplit;
		sl=c.sl;
		fileFormat=c.fileFormat;
		mzXMLFilter=c.mzXMLFilter;
		distArea=c.distArea;
		strcpy(formula,c.formula);

		centroid = c.centroid;
		staticSN = c.staticSN;
    xml = c.xml;
		reducedOutput = c.reducedOutput;

    //ppMatch=c.ppMatch;
    //ppWin=c.ppWin;
    //noiseMatch=c.noiseMatch;
    //noiseWindow=c.noiseWindow;
    ppm=c.ppm;
    //sna=c.sna;
		boxcar=c.boxcar;
		boxcarFilter=c.boxcarFilter;

    //rawAvg=c.rawAvg;
    //rawAvgWidth=c.rawAvgWidth;
    //rawAvgCutoff=c.rawAvgCutoff;

    strcpy(rawFilter,c.rawFilter);
  }
  return *this;
}

void CHardklorSetting::clearVariant(){
  delete variant;
  variant = new vector<CHardklorVariant>;
}

void CHardklorSetting::out(char *s){
  sprintf(s,"minCh:%d maxCh:%d d:%d p:%d s:%d corr:%lf sn:%lf res:%d,%lf win:%lf sl:%d v:%d\n",minCharge,maxCharge,depth,peptide,smooth,corr,sn,msType,res400,winSize,sl,variant->size());
}
