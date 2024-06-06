#include "CnpxQuanticSummary.h"

using namespace std;

CnpxQuanticSummary::CnpxQuanticSummary() {
  version.clear();
  options.clear();
}

void CnpxQuanticSummary::write(FILE* f){
  size_t i;

  fprintf(f, "<quantic_summary version=\"%s\"", version.c_str());
  fprintf(f, " options=\"%s\"", options.c_str());
  fprintf(f, ">\n");

  for (i = 0; i<inputfile.size(); i++) inputfile[i].write(f);

  fprintf(f, "</interprophet_summary>\n");
}