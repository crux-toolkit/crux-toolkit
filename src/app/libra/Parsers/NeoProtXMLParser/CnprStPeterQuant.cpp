#include "CnprStPeterQuant.h"

using namespace std;

CnprStPeterQuant::CnprStPeterQuant(){
  dSI=0;
  dSIn=0;
  SI=0;
  SIn=0;
  dCounts=0;
  counts=0;
  dNSAF=0;
  NSAF=0;
  ng=0;
  ngC=0;
}

void CnprStPeterQuant::write(FILE* f, int tabs){
  //required
  string el = "StPeterQuant";
  if (dSIn==0 && SIn==0) NPRerrMsg(el, "dSIn or SIn");

  NPRprintTabs(f, tabs);
  fprintf(f, "<StPeterQuant");
  if(dSI!=0)fprintf(f, " dSI=\"%.10lf\"", dSI);
  if(dSIn!=0) fprintf(f, " dSIn=\"%.10lf\"", dSIn);
  if (SI != 0)fprintf(f, " SI=\"%.10lf\"", SI);
  if (SIn != 0) fprintf(f, " SIn=\"%.10lf\"", SIn);
  if (dCounts != 0)fprintf(f, " dCounts=\"%.4lf\"", dCounts);
  if (counts != 0) fprintf(f, " counts=\"%.4lf\"", counts);
  if (dNSAF != 0)fprintf(f, " dNSAF=\"%.10lf\"", dNSAF);
  if (NSAF != 0) fprintf(f, " NSAF=\"%.10lf\"", NSAF);
  if (ng != 0)fprintf(f, " ng=\"%.10lf\"", ng);
  if (ngC != 0) fprintf(f, " ngC=\"%.10lf\"", ngC);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<StPeterQuant_peptide.size(); j++) StPeterQuant_peptide[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</StPeterQuant>\n");

}
