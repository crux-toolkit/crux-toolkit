#ifndef _CSPLITSPECTRUM_H
#define _CSPLITSPECTRUM_H

#include "CHardklorSetting.h"
#include "CSpecAnalyze.h"
#include "Spectrum.h"
#include <vector>

using namespace std;

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

	vector<CSpecAnalyze> *setA;
	vector<CSpecAnalyze> *setB;
	vector<CSpecAnalyze> *finalAnalysis;

	vector<float> *s2n;

	vector<int> *aIndex;
	vector<int> *bIndex;

};

#endif

