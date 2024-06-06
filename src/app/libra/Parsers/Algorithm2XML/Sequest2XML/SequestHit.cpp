#include "SequestHit.h"

SequestHit::SequestHit() {
  dMass = 0.0;
  dXC = 0.0;
  dDeltCn = 0.0;
  dSp = 0.0;
  cAA1 = '\0';
  cAA2 = '\0';
  szProt[0] = '\0';
  szPlainPep[0] = '\0';
  szSubPep[0] = '\0';
  szDSite[0] = '\0';
  szDup[0] = '\0';
  iRankSp = 0;
  iIon = 0;
  iTot = 0;
  iDeltCnIdxDiff = -1;
  dSpecialDeltCn = 0;
}

void SequestHit::writeOutFile(FILE* outFile, double deltaCn) {
  fprintf(outFile, " %d.   %2d /%3d  %4.4f  %1.4f   %1.4f  %3.1f  %2d/%2d  %s %s %c.%s.%c\n", 
	  iRank, iRank, iRankSp, dMass, deltaCn, dXC, dSp,
	  iIon, iTot, szProt, szDup, cAA1, szSubPep, cAA2);

}
