#include "CnpxInteractSummary.h"

using namespace std;

void CnpxInteractSummary::write(FILE* f){
  size_t i;

  fprintf(f, "<interact_summary filename=\"%s\"", filename.c_str());
  fprintf(f, " directory=\"%s\"", directory.c_str());
  fprintf(f, ">\n");

  for (i = 0; i<inputfile.size(); i++) inputfile[i].write(f);

  fprintf(f, "</interact_summary>\n");
}