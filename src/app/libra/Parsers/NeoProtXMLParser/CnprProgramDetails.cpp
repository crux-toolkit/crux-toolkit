#include "CnprProgramDetails.h"

using namespace std;

void CnprProgramDetails::write(FILE* f, int tabs){
  //required
  string el = "program_details";
  if (analysis.empty()) NPRerrMsg(el, "analysis");
  if (time.date.year==0) NPRerrMsg(el, "time");

  NPRprintTabs(f, tabs);
  fprintf(f, "<program_details");
  fprintf(f, " analysis=\"%s\"", analysis.c_str());
  fprintf(f, " time=\"%s\"", time.write().c_str());
  if (!version.empty()) fprintf(f, " type=\"%s\"", version.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<proteinprophet_details.size(); j++) proteinprophet_details[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</program_details>\n");

}
