#ifndef _CNOISEREDUCTION_H
#define _CNOISEREDUCTION_H

#include "MSReader.h"
#include "Spectrum.h"
#include "CHardklorSetting.h"
#include <cmath>
#include <iostream>
#include <deque>

#define GC 5.5451774444795623

class CNoiseReduction {

public:

  //Constructors and Destructors
  CNoiseReduction();
  CNoiseReduction(MSReader* msr, CHardklorSetting& hs);
  ~CNoiseReduction();

  //Functions
  double CParam(Spectrum& sp, int tot=1);
  double calcFWHM(double mz);
  void FirstDerivativePeaks(Spectrum& sp, int winSize);
  bool DeNoise(Spectrum& sp);
  //bool DeNoise(Spectrum& sp, vector<Spectrum>& vs, int pivot, bool findPeaks=false);
  //bool DeNoise(Spectrum& sp, deque<Spectrum>& vs, int pivot, bool findPeaks=false);
  bool DeNoiseB(Spectrum& sp);
  bool DeNoiseC(Spectrum& sp);
  bool DeNoiseD(Spectrum& sp);
  int NearestPeak(Spectrum& sp, double mz);
  bool ScanAverage(Spectrum& sp, char* file, int width, float cutoff);
  bool NewScanAverage(Spectrum& sp, char* file, int width, float cutoff, int scanNum=0);
  //bool ScanAverage(Spectrum& sp, vector<Spectrum>& vs, int pivot, float cutoff, double cp=0.0);
  //bool ScanAverage(Spectrum& sp, deque<Spectrum>& vs, int pivot, float cutoff);
  //bool ScanAverageBuffered(Spectrum& sp, char* file, int width, float cutoff, int scanNum=0);
  bool ScanAveragePlusDeNoise(Spectrum& sp, char* file, int width, float cutoff, int scanNum=0);
  bool NewScanAveragePlusDeNoise(Spectrum& sp, char* file, int width, float cutoff, int scanNum=0);

  int pos;

private:
  //Functions
  
  
  //Data Members
  //int pos;
  int posA;
  char lastFile[256];
  CHardklorSetting cs;
  MSReader* r;
  deque<Spectrum> s;
  deque<Spectrum> bs;

	/*
	  __int64 startTime;
    __int64 stopTime;
    __int64 loadTime;
    __int64 analysisTime;
    __int64 timerFrequency;
    __int64 tmpTime1;
    __int64 tmpTime2;
    #define getExactTime(a) QueryPerformanceCounter((LARGE_INTEGER*)&a)
    #define getTimerFrequency(a) QueryPerformanceFrequency((LARGE_INTEGER*)&a)
    #define toMicroSec(a) (a)
    #define timeToSec(a,b) (a/b)
		*/

};

#endif

