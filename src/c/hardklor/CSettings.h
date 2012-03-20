#ifndef _CSETTINGS_H
#define _CSETTINGS_H

#include <vector>

using namespace std;

enum ScanType{
  Zoom,
  UltraZoom,
  IonSpec2,
  Other
};

typedef struct {
  int code;
  double value;
} varType;

typedef struct {
  char infile[256];
  char outfile[256];
  ScanType st;
} FileName;

class CSettings {
 public:
  bool QAR;
  bool ROC;
  int chlor;
  int selen;
  int smooth;
  int maxPep;
  int maxCharge;
  int scan;
  int maxDepth;
  double S2N;
  double signal;
  double corrThresh;
  double windowLower;
  double windowUpper;
  vector<FileName*> files;
  vector<varType> variants;
  CSettings();
  CSettings(char*);
  char enrichAtom[3];
  double enrich;
  int enrichTimes;
 protected:
 private:
  void readFile(char*);
};



#endif
