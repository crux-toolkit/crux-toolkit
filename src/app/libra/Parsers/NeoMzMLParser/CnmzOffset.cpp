#include "CnmzOffset.h"

using namespace std;

void CnmzOffset::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "offset";
  if (idRef.empty()) NMZerrMsg(el, "idRef");

  NMZprintTabs(f, tabs);
  fprintf(f, "<offset");
  fprintf(f, " idRef=\"%s\"", idRef.c_str());
  if (scanTime>0) fprintf(f, " scanTime=\"%.4lf\"", scanTime);
  if (!spotID.empty()) fprintf(f, " spotID=\"%s\"", spotID.c_str());
  fprintf(f,">");
  fprintf(f,"%s",content.c_str());
  fprintf(f, "</offset>\n");



}

