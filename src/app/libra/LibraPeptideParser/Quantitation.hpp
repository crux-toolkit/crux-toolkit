#ifndef QUANTITATION_HPP
#define QUANTITATION_HPP

#include "LibraConditionHandler.hpp"

#include <fstream>	// file output
#include <iomanip>	// setw
#include <iostream>
#include <vector>
#include <math.h>	// fabs
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Parsers/mzParser/cached_ramp.h" // deals with inefficient use of RAMP file open etc
#include "typedefs.hpp"
#include <map>


using std::vector;
using std::pair;
using std::make_pair;
using std::map;

#define CENTROID 0.5

class Quantitation {
private:
  LibraConditionHandler* pLibraConditionHandler;

  std::ofstream m_fout;
  double	m_RT;
  char* mzXMLFile;
  mzParser::RAMPFILE  *pFI;

  bool writeToOutFile;
  bool printToStdOut;

  int startScanNum, endScanNum;

  vff	m_maxima;		// mz-Int in the reporter ions mz interval
  vff	m_maxima_absolute;      // same as m_maxima, but never normalized
  vff	m_centroidPeaks;
  vff	m_profilePeaks;
  ff	m_massInt;
  vff::iterator	pos, pos2;

  // making maps to hold:
  // key = scan number, value = vector of target masses
  map<int, vector<double> > mapScanNumAndTargetMasses;

  // key = scan number, value = vector of absolute intensities
  map<int, vector<double> > mapScanNumAndAbsoluteIntensities;

  // key = scan number, value = vector of normalized intensities
  map<int, vector<double> > mapScanNumAndNormalizedIntensities;

  int   m_centroidingIteration; // number of times the centroiding process

  double m_tic;

public:
  Quantitation(LibraConditionHandler*, char*, mzParser::RAMPFILE*);

  ~Quantitation();

  int findMaxima();
  int centroid(int);
  int centroidWeighted(int);
  int isotopeCorrection();
  int plotCentroid(int);
  int printMaxima();
  int normalizeRatio();
  int printRatio(int);
  int reset();
  void setRT(double);
  vff getMaxima();
  int getCentroidingIteration();
  int getIsToNormalize();
  int calculateIntensities();
  int calculateIntensities(int startScan, int endScan);
  int calculateIntensities(mzParser::RAMPFILE *, mzParser::ramp_fileoffset_t *, int, int, int);

  int openExistingOutFile(char*);
  void printVFFToStdout(vff);

  vector<double> getTargetMasses(int scanNumber);
  vector<double> getAbsoluteIntensities(int scanNumber);
  vector<double> getNormalizedIntensities(int scanNumber);

  int getNumberOfChannels(int scanNumber);
  int getStartScanNumber();
  int getEndScanNumber();

  void setWriteToOutFile(bool);
  void setPrintToStdOut(bool); 

};

#endif
