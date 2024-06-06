#ifndef _CNPXUIPSM_H
#define _CNPXUIPSM_H

#include "CnpxSpectrumQuery.h"
#include <string>
#include <vector>

typedef struct npxUIMod{
  int pos;
  double mass;
  double diffMass;
} npxUIMod;

typedef struct npxUIModification{
  std::string modifiedPeptide;
  std::vector<npxUIMod> mods;
} npxUIModification;

typedef struct npxUIProtein{
  char nextAA;
  char prevAA;
  std::string protein;
} npxUIProtein;

typedef struct npxUIScore{
  std::string name;
  std::string value;
} npxUIScore;

typedef struct npxUIProbability{
  double probability;
  std::vector<npxUIScore> parameters;
} npxUIProbability;

class CnpxUIPSM {
public:
  int assumed_charge;
  double calcNeutralMass;
  bool hasPeptideProphet;
  bool hasIProphet;
  npxUIProbability iProphet;
  npxUIModification mod;
  bool modified;
  std::string peptide;
  npxUIProbability peptideProphet;
  double precursorNeutralMass;
  std::vector<npxUIProtein> proteins;
  double retentionTime;
  int scanNumber;
  std::vector<npxUIScore> scores;
  std::string spectrumID;
  int xlType;

  void setPSM(CnpxSpectrumQuery& s);

private:

};

#endif 