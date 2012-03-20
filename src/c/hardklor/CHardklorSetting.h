#ifndef _CHARDKLORSETTING_H
#define _CHARDKLORSETTING_H

#include <string>
#include <vector>

#include "HardklorTypes.h"
#include "CHardklorVariant.h"

using namespace std;

/* 
   Defined as a class instead of a struct so
   that default parameters are used. Data access
   will be identical to that of a struct.
*/
class CHardklorSetting {

 public:
  //Constructors & Destructors:
  CHardklorSetting();
  CHardklorSetting(const CHardklorSetting&);
  ~CHardklorSetting();

  //Methods:
  CHardklorSetting& operator=(const CHardklorSetting&);
  void clearVariant();
  void out(char *s);

  //Data Mebers:
  bool centroid;    //spectrum data is centroided
  bool distArea;    //report distribution area instead of base peak intensity
  bool iAnalysis;   //intersect analysis(true) or union analysis(false)
  bool noBase;      //No base molecule - perform analysis with only averagine variant models
  bool noSplit;     //analyze entire spectrum at once
  //bool rawAvg;      //use averaged raw scans
  bool skipZero;    //ignore zero intensity data points
  bool staticSN;    //for sna=THRASH; assume one noise level for entire spectrum
  bool xml;         //output is in xml
  
  
  int depth;        //maximum number of overlapping peptides
  int maxCharge;    //max charge state to search for
  int minCharge;    //min charge state to search for
  int peptide;			//maximum peptide models to analyze at a single time
  int ppMatch;      //pre-processing matches. m/z must be observed this amount across ppWin
  int ppWin;        //pre-processing window size (1 = +/-1 scan)
  //int noiseMatch;   //for sna=PP; number of matches required for real peaks
  //int noiseWindow;  //for sna=PP; Size of window over which scans are analyzed
  //int rawAvgCutoff; //Noise cutoff intensity for averaged raw scans
  //int rawAvgWidth;  //Number of scans on either side of target to average (1 = +/-1 scan)
  int sl;           //sensitivity level
  int smooth;       //Savitsky-Golay smoothing window size
  int sna;          //Signal-to-noise algorithm; 0=THRASH, 1=Persistent peaks (PP)

  double corr;      //correlation threshold
  double ppm;       //ppm tolerance of m/z values to match across scans
  double res400;    //resolution at m/z 400
  double sn;        //for sna=THRASH; signal-to-noise ratio threshold
  double snWindow;  //for sna=THRASH; bin size over which local noise is computed
  double winSize;   //maximum window size for analysis

  sInt scan;        //scan range to analyze
  sDouble window;   //m/z range to analyze

  char chargeMode;        //charge determination function to use
  char formula[64];       //non-averagine model formula
  char inFile[256];       //input file name
  char outFile[256];      //output file name
  char rawFilter[256];    //Filter which spectra from raw files are analyzed
  char MercuryFile[256];  //mercury data file to use
  char HardklorFile[256]; //hardklor data file to use

  specType msType;                    //Type of mass spectrometer used to acquire data
  hkAlgorithm algorithm;              //Deconvolving algorithm to use
  vector<CHardklorVariant> *variant;  //Variants to make to averagine
  MSFileFormat fileFormat;            //File format
  MSSpectrumType mzXMLFilter;         //Filter for mzXML files

 protected:
 private:

};

#endif
