#ifndef _HARDKLORTYPES_H
#define _HARDKLORTYPES_H

#include "MSToolkitTypes.h"

#include <string>
#include <vector>

using namespace std;
using namespace MSToolkit;

enum specType {
	OrbiTrap,
	TOF,
	QIT,
	FTICR
};

enum hkAlgorithm {
	Basic,
	SemiComplete,
	SemiCompleteFast,
	Dynamic,
	DynamicSemiComplete,
	SemiSubtractive,
	FewestPeptides,
	FewestPeptidesChoice,
	FastFewestPeptides,
	FastFewestPeptidesChoice
};


typedef struct {
  char atom[3];
  int isotope;
  double ape;
} sEnrich;

typedef struct {
  string molecule;
  int iLower;
  int iUpper;
} sMolecule;

typedef struct {
  int iLower;
  int iUpper;
} sInt;

typedef struct{
  double dLower;
  double dUpper;
} sDouble;

typedef struct{
  float fLower;
  float fUpper;
} sFloat;

typedef struct{
  int iValue;
  double dValue;
} sID;

typedef struct{
  int atomNum;
  int isotope;
  double ape;
} sEnrichMercury;

typedef struct{
  double mz;
  float intensity;
  int index;
} sSplit;

//Do this better
enum ScanType{
  Zoom,
  UltraZoom,
  IonSpec2,
  Other
};

typedef struct {
  char id[5];
  double mz;
  double monoMass;
  double shft;
  double abun;
  int charge;
  char seq[31];
  int C;
  int H;
  int O;
  int N;
  int S;
  vector<sID> *enrich;
} peps;

typedef struct pepHit{
	int basePeakIndex;
	int charge;
	int lowIndex;
	int highIndex;
	int variantIndex;
	float intensity;
	float area;
	double massShift;
	double lowMZ;
	double highMZ;
	double monoMass;
	double corr;
	char averagine[32];
} pepHit;

typedef struct mercuryModel{
	float area;
	int size;
	double zeroMass;
	Peak_T* peaks;
} mercuryModel;

#endif
