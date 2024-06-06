#include "CnprParameter.h"

using namespace std;

void CnprParameter::write(FILE* f, int tabs){
  //required
  string el = "parameter";
  if (name.empty()) NPRerrMsg(el, "name");
  if (value.empty()) NPRerrMsg(el, "value");

  NPRprintTabs(f, tabs);
  fprintf(f, "<parameter");
  fprintf(f, " name=\"%s\"", name.c_str());
  fprintf(f, " value=\"%s\"", value.c_str());
  if (!type.empty()) fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f, "/>\n");

}
