#include "CnmzCvParam.h"

using namespace std;

void CnmzCvParam::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "cvParam";
  if (accession.empty()) NMZerrMsg(el, "accession");
  if (cvRef.empty()) NMZerrMsg(el, "cvRef");
  if (name.empty()) NMZerrMsg(el, "name");

  NMZprintTabs(f, tabs);
  fprintf(f, "<cvParam");
  fprintf(f, " accession=\"%s\"", accession.c_str());
  fprintf(f, " cvRef=\"%s\"", cvRef.c_str());
  fprintf(f, " name=\"%s\"", name.c_str());
  if (!value.empty())fprintf(f, " value=\"%s\"", value.c_str());
  if (!unitAccession.empty()) fprintf(f, " unitAccession=\"%s\"", unitAccession.c_str());
  if (!unitCvRef.empty())fprintf(f, " unitCvRef=\"%s\"", unitCvRef.c_str());
  if (!unitName.empty())fprintf(f, " unitName=\"%s\"", unitName.c_str());
  
  fprintf(f, "/>\n");



}

