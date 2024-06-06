#include "CnmzSoftwareRef.h"

using namespace std;

void CnmzSoftwareRef::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "softwareRef";

  NMZprintTabs(f, tabs);
  fprintf(f, "<softwareRef");
  fprintf(f, " ref=\"%s\"", ref.c_str());
  fprintf(f, "/>\n");

}

