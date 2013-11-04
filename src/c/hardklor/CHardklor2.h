#ifndef _CHARDKLOR2_H
#define _CHARDKLOR2_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <cmath>

#include "MSObject.h"
#include "MSReader.h"
#include "Spectrum.h"
#include "CAveragine.h"
#include "CMercury8.h"
#include "CHardklor.h"
#include "CModelLibrary.h"

#ifdef _MSC_VER

#else
#include <sys/time.h>
#endif

using namespace std;

class CHardklor2{

 public:
  //Constructors & Destructors:
	CHardklor2(CAveragine *a, CMercury8 *m, CModelLibrary *lib);
  ~CHardklor2();

  //Methods:
	void Echo(bool b);
  int GoHardklor(CHardklorSetting sett);

 protected:

 private:
  //Methods:
	int			BinarySearch(Spectrum& s, double mz, bool floor);
	double	CalcFWHM(double mz,double res,int iType);
	void		Centroid(Spectrum& s, Spectrum& out);
	bool		CheckForPeak(vector<Result>& vMR, Spectrum& s, int index);
  int			CompareData(const void*, const void*);
  double	LinReg(vector<float>& mer, vector<float>& obs);
	bool		MatchSubSpectrum(Spectrum& s, int peakIndex, pepHit& pep);
	double	PeakMatcher(vector<Result>& vMR, Spectrum& s, double lower, double upper, double deltaM, int matchIndex, int& matchCount, int& indexOverlap, vector<int>& vMatchIndex, vector<float>& vMatchIntensity);
	double	PeakMatcherB(vector<Result>& vMR, Spectrum& s, double lower, double upper, double deltaM, int matchIndex, int& matchCount, vector<int>& vMatchIndex, vector<float>& vMatchIntensity);
	void		QuickCharge(Spectrum& s, int index, vector<int>& v);
	void		QuickHardklor(Spectrum& s, vector<pepHit>& vPeps);
	void		RefineHits(vector<pepHit>& vPeps, Spectrum& s);
  void		WritePepLine(pepHit& ph, Spectrum& s, FILE* fptr, int format=0); 
  void		WriteScanLine(Spectrum& s, FILE* fptr, int format=0); 

	static int CompareBPI(const void *p1, const void *p2);

  //Data Members:
	CHardklorSetting cs;
	CAveragine *averagine;
	CMercury8 *mercury;
	CModelLibrary* models;
	CPeriodicTable* PT;
	Spectrum mask;
	bool bEcho;

  //Temporary Data Members:
  char bestCh[200];
  double BestCorr;
  int CorrMatches;
  int ExtraPre;
  int ExtraObs;

	int SSIterations;

  //For accurate timing of Hardklor
  #ifdef _MSC_VER
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
  #else
    timeval startTime;
    timeval stopTime;
    uint64_t loadTime;
    uint64_t splitTime;
    uint64_t analysisTime;
    uint64_t tmpTime1;
    uint64_t tmpTime2;
    int timerFrequency;
    #define getExactTime(a) gettimeofday(&a,NULL)
    #define toMicroSec(a) a.tv_sec*1000000+a.tv_usec
    #define getTimerFrequency(a) (a)=1
    #define timeToSec(a,b) (a/1000000)
  #endif

};

#endif
