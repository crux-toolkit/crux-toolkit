#include "MSReader.h"
#include "Spectrum.h"
#include "CHardklorSetting.h"
#include <cmath>
#include <iostream>
#include <deque>

#define GC 5.5451774444795623


using namespace MSToolkit;

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
  bool NewScanAverage(Spectrum& sp, char* file, int width, /*float cutoff,*/ int scanNum=0);
  //bool ScanAverage(Spectrum& sp, vector<Spectrum>& vs, int pivot, float cutoff, double cp=0.0);
  //bool ScanAverage(Spectrum& sp, deque<Spectrum>& vs, int pivot, float cutoff);
  //bool ScanAverageBuffered(Spectrum& sp, char* file, int width, float cutoff, int scanNum=0);
  bool ScanAveragePlusDeNoise(Spectrum& sp, char* file, int width, /*float cutoff,*/ int scanNum=0);
  bool NewScanAveragePlusDeNoise(Spectrum& sp, char* file, int width, /*float cutoff,*/ int scanNum=0);

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

};

