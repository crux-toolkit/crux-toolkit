#include "CnpxParameter.h"

using namespace std;

void CnpxParameter::write(FILE* f, int tabs) {
  string el = "parameter";
  if (name.empty()) NPXerrMsg(el, "name");
  if (value.empty()) NPXerrMsg(el, "value");

  int t = tabs;
  if (t>-1) t++;

  NPXprintTabs(f, tabs);
  fprintf(f, "<parameter name=\"%s\" value=\"%s\"", name.c_str(), value.c_str());
  if (!type.empty()) fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f, "/>\n");

}