#include "CnmzFileChecksum.h"

using namespace std;

void CnmzFileChecksum::write(FILE* f, int tabs, bool iterative){

  NMZprintTabs(f, tabs);
  fprintf(f, "<fileChecksum>");
  fprintf(f, "%s", content.c_str());
  fprintf(f, "</fileChecksum>\n");

}

