#ifndef _CSPLITSPECTRUM_H
#define _CSPLITSPECTRUM_H

#include "CHardklorSetting.h"
#include "CSpecAnalyze.h"
#include "MSToolkit/Spectrum.h"
#include <vector>

class CSplitSpectrum {
 public:
  //Constructors & Destructors
  CSplitSpectrum(Spectrum *spec, CHardklorSetting& sett);
	~CSplitSpectrum();
  
  //Fuctions
  int getNumWindows();
	CSpecAnalyze getWindow(int index);
	void IntersectionAnalysis();
	void MakeAnalysis(double winSize);
	void NoSplitAnalysis();
	void OverlappingAnalysis(double gapSize);
	void SinglePassAnalysis(double gapSize);
	void UnionAnalysis();

	void SetAveragine(CAveragine *a);
	void SetMercury(CMercury8 *m);

	void Centroid(Spectrum& s);

	void NewSNPass(double gapSize);


 private:
  //Data members
  Spectrum *wholeSpec;
	CHardklorSetting userParams;

	CSpecAnalyze goodPeaks;
	CAveragine *averagine;
	CMercury8 *mercury;

  std::vector<CSpecAnalyze> *setA;
  std::vector<CSpecAnalyze> *setB;
  std::vector<CSpecAnalyze> *finalAnalysis;

  std::vector<float> *s2n;

  std::vector<int> *aIndex;
  std::vector<int> *bIndex;

};

#endif

