#include "CnmzInstrumentConfiguration.h"

using namespace std;

void CnmzInstrumentConfiguration::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "instrumentConfiguration";
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<instrumentConfiguration");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " scanSettingsRef=\"%s\"", scanSettingsRef.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);
  for (size_t a = 0; a<componentList.size(); a++) componentList[a].write(f, t, iterative);
  for (size_t a = 0; a<softwareRef.size(); a++) softwareRef[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</instrumentConfiguration>\n");

}

