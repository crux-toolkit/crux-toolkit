#ifndef _SEQUESTHIT_H_
#define _SEQUESTHIT_H_
#include "Common/constants.h"
#include <stdio.h>
class SequestHit {
  friend class SequestOut;
  friend class Out2XML;
  friend class CombineOut;
 public:
  SequestHit();
  ~SequestHit() {};
  void writeOutFile(FILE* outFile, double deltaCn);

 private:
  double dMass;
  double dXC;
  double dDeltCn;
  double dSp;
  char cAA1;
  char cAA2;
  char szProt[SIZE_PEP];
  char szPlainPep[SIZE_PEP];
  char szSubPep[SIZE_PEP];
  char szDSite[SIZE_PEP];
  char szDup[SIZE_PEP];
  int  iRankSp;
  int  iRank;
  int  iIon;
  int  iTot;
  double  dSpecialDeltCn;
  int  iDeltCnIdxDiff;
};

#endif //_SEQUESTHIT_H_
