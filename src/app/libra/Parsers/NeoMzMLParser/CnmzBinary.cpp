#include "CnmzBinary.h"

using namespace std;

void CnmzBinary::write(FILE* f, int tabs, bool iterative){
  NMZprintTabs(f, tabs);
  fprintf(f, "<binary>");
  fprintf(f, "%s",content.c_str());
  fprintf(f, "</binary>\n");
}

