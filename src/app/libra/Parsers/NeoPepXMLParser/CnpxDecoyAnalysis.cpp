#include "CnpxDecoyAnalysis.h"

using namespace std;

void CnpxDecoyAnalysis::write(FILE* f){
  size_t i;

  fprintf(f, "<decoy_analysis>\n");
  for (i = 0; i<point.size(); i++) point[i].write(f);
  fprintf(f, "</decoy_analysis>\n");
}