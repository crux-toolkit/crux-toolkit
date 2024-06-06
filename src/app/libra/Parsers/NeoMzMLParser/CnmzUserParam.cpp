#include "CnmzUserParam.h"

using namespace std;

void CnmzUserParam::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "userParam";
  if (name.empty()) NMZerrMsg(el, "name");

  NMZprintTabs(f, tabs);
  fprintf(f, "<userParam");
  fprintf(f, " name=\"%s\"", name.c_str());
  if (!type.empty())fprintf(f, " type=\"%s\"", type.c_str());
  if (!value.empty())fprintf(f, " value=\"%s\"", value.c_str());
  if (!unitAccession.empty()) fprintf(f, " unitAccession=\"%s\"", unitAccession.c_str());
  if (!unitCvRef.empty())fprintf(f, " unitCvRef=\"%s\"", unitCvRef.c_str());
  if (!unitName.empty())fprintf(f, " unitName=\"%s\"", unitName.c_str());

  fprintf(f, "/>\n");



}

