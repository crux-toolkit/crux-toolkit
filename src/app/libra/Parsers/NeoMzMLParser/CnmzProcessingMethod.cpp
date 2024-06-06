#include "CnmzProcessingMethod.h"

using namespace std;

void CnmzProcessingMethod::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "processingMethod";
  if (order<-999) NMZerrMsg(el, "order");
  if (softwareRef.empty()) NMZerrMsg(el, "order");

  NMZprintTabs(f, tabs);
  fprintf(f, "<processingMethod");
  fprintf(f, " order=\"%d\"", order);
  fprintf(f, " softwareRef=\"%s\"", softwareRef.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</processingMethod>\n");

}

