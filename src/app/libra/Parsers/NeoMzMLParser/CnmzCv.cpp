#include "CnmzCv.h"

using namespace std;

void CnmzCv::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "cv";
  if (URI.empty()) NMZerrMsg(el, "URI");
  if (fullName.empty()) NMZerrMsg(el, "fullName");
  if (id.empty()) NMZerrMsg(el, "id");

  NMZprintTabs(f, tabs);
  fprintf(f, "<cv");
  fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, " fullName=\"%s\"", fullName.c_str());
  if (!version.empty())fprintf(f, " version=\"%s\"", version.c_str());
  fprintf(f, " URI=\"%s\"", URI.c_str());
  fprintf(f, "/>\n");

}

