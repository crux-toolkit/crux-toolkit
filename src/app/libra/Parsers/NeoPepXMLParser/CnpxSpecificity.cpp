#include "CnpxSpecificity.h"

using namespace std;

CnpxSpecificity::CnpxSpecificity(){
  cut.clear();
  min_spacing=0;
  no_cut.clear();
  sense.clear();
}

void CnpxSpecificity::write(FILE* f, int tabs){

  string el = "specificity";
  if (sense.empty()) NPXerrMsg(el, "sense");
  if (cut.empty()) NPXerrMsg(el, "cut");

  NPXprintTabs(f, tabs);
  fprintf(f, "<specificity cut=\"%s\"", cut.c_str());
  if (!no_cut.empty()) fprintf(f, " no_cut=\"%s\"", no_cut.c_str());
  fprintf(f, " sense=\"%s\"", sense.c_str());
  if (min_spacing>0) fprintf(f, " min_spacing=\"%u\"",min_spacing);
  fprintf(f, "/>\n");

}