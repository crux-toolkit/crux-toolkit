#include "CnmzIndexListOffset.h"

using namespace std;

void CnmzIndexListOffset::write(FILE* f, int tabs, bool iterative){

  NMZprintTabs(f, tabs);
  fprintf(f, "<indexListOffset>");
  fprintf(f, "%s", content.c_str());
  fprintf(f, "</indexListOffset>\n");

}

