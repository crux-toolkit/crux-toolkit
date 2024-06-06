#include "CnmzIndexedmzML.h"

using namespace std;

void CnmzIndexedmzML::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "CnmzIndexedmzML";

  NMZprintTabs(f, tabs);
  fprintf(f, "<indexedmzML");
  fprintf(f, " xmlns=\"%s\"", xmlns.c_str());
  fprintf(f, " xmlns:xsi=\"%s\"", xmlns_xsi.c_str());
  fprintf(f, " xsi:schemaLocation=\"%s\"", xsi_schemaLocation.c_str());
  fprintf(f, ">\n");

  int t=tabs;
  if(t>-1) t++;
  mzML->write(f,t,iterative);

  if (iterative) return;

  indexList.write(f,t,iterative);
  indexListOffset.write(f,t,iterative);
  fileChecksum.write(f,t,iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</indexedmzML>\n");

}