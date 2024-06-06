#include "CnprStPeterQuantPeptide.h"

using namespace std;

CnprStPeterQuantPeptide::CnprStPeterQuantPeptide(){
  charge=0;
  dSI=0;
  SI=0;
  dSC=0;
  SC=0;
}

void CnprStPeterQuantPeptide::write(FILE* f, int tabs){
  //required
  string el = "StPeterQuant_peptide";
  if (sequence.empty()) NPRerrMsg(el, "sequence");
  if (charge == 0)  NPRerrMsg(el, "charge");

  NPRprintTabs(f, tabs);
  fprintf(f, "<StPeterQuant_peptide");
  fprintf(f, " sequence=\"%s\"", sequence.c_str());
  fprintf(f, " charge=\"%d\"", charge);
  if (dSI != 0) fprintf(f, " dSI=\"%.1lf\"", dSI);
  if (SI != 0) fprintf(f, " SI=\"%.1lf\"", SI);
  if (dSC != 0) fprintf(f, " dSC=\"%.4lf\"", dSC);
  if (SC != 0) fprintf(f, " SC=\"%.4lf\"", SC);
  fprintf(f, "/>\n");

}
