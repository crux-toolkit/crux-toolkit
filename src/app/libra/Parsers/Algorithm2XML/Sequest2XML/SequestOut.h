#ifndef _SEQUESTOUT_H_
#define _SEQUESTOUT_H_
#include "Common/constants.h"
#include "Common/Array.h"
#include "SequestHit.h"
#include "SequestParams.h"
#include "Common/TPPVersion.h"

class SequestOut {
  friend class SequestHit;
  friend class Out2XML;
  friend class CombineOut;
 public:
  SequestOut();
  SequestOut(const SequestOut& out);
  ~SequestOut();
  int getNumHits();
  void insertNextHit(SequestHit* hit);
  void getDeltaCn();
  SequestHit* getHitByIndex(int idx);
  void writeOutFile(FILE* outFile, SequestParams* params);

 private:
  char szFileName[SIZE_FILE];
  char szBaseFileName[SIZE_FILE];

  char szMod[SIZE_FILE];
  char szDatabase[SIZE_FILE];
  double dAMass;
  double dMass;

  double dMass1;
  double dMass2;
  double dMass3;
  double dMass4;
  double dMass5;
  double dMass6;
  double dMassCT;
  double dMassNT;

  int  iMassType;
  int  bNucDb;
  Array<SequestHit*> sequestHits_;
};

#endif //_SEQUESTOUT_H_
