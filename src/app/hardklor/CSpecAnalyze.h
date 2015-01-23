#ifndef _CSPECANALYZE_H
#define _CSPECANALYZE_H

#include "CAveragine.h"
#include "CMercury8.h"
#include "HardklorTypes.h"
#include "SpecAnalyzeSupport.h"
#include "Spectrum.h"
#include "S2N.h"
#include "CHardklorSetting.h"
#include "FFT-HK.h"
#include <vector>
#include <cmath>
using namespace std;

//the following is 8*ln(2)
#define GAUSSCONST 5.5451774444795623

class CSpecAnalyze {
 public:
  //Constructors & Destructors
  CSpecAnalyze();
  CSpecAnalyze(const CSpecAnalyze&);
	~CSpecAnalyze();

  //Overloaded operators
  CSpecAnalyze& operator=(const CSpecAnalyze&);

  //Functions
  void BuildMismatchArrays();
	void chargeState();
  void clear();
	//void FindCharge();
  int  FindPeaks();
	int  FindPeaks(Spectrum& s, int start, int stop);
  //void FindPeptides();
  void MakePredictions(vector<CHardklorVariant>& var);
  int  PredictPeptides();
  void removePeaksBelowSN();
  void setSpectrum(const Spectrum& s);
	void setAveragine(CAveragine *a);
	void setMercury(CMercury8 *m);
	void setParams(const CHardklorSetting& sett);
	void TraditionalCharges();

  //Data Members
  float                      basePeak;     //most intense peak in peaks spectrum; used for SN cutoff
  vector<int>                *charges;
  bool                       manyPeps;
  int                        mismatchSize;
  Spectrum                   peaks;        //peaks found in the spectrum
  Spectrum                   peptide;      //peaks that are potential peptides
  vector<CPeakPrediction>    *predPeak;
	vector<CPeptidePrediction> *predPep;
  float                      S2NCutoff;	

 protected:
 private:

  //Functions
	double calcFWHM(double mz);
  void   FirstDerivativePeaks(Spectrum&,int);
	void   FirstDerivativePeaks(Spectrum&,int,int,int);
	double InterpolateMZ(Peak_T& p1, Peak_T& p2, double halfMax);

  //Data members
  CAveragine       *averagine;
	CMercury8        *mercury;
  Spectrum         spec;
	CHardklorSetting userParams;
	
};


#endif

