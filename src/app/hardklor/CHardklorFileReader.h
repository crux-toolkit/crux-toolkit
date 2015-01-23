#include "MSReader.h"
#include "Spectrum.h"
#include "CHardklorSetting.h"
#include "CNoiseReduction.h"
#include <cmath>
#include <iostream>
#include <deque>

#define GC 5.5451774444795623

typedef struct cpar{
  int scan;
  double c;
} cpar;

class CHardklorFileReader {
public:
  CHardklorFileReader();
  CHardklorFileReader(CHardklorSetting& hs);
  ~CHardklorFileReader();

  int getPercent();
  //bool getProcessedSpectrum(Spectrum& s);
  //bool getProcessedSpectrumB(Spectrum& s);
  //bool initBuffers();
  //bool getNRSpectrum(Spectrum& s, int pos);

private:

  //Data Members
  int posNative;
  int posNR;
  CHardklorSetting cs;
  CNoiseReduction nr;
  MSReader r;
  deque<Spectrum> sNative;
  deque<Spectrum> sNR;

  deque<cpar> dc;

  vector<int> vScans;
  
};