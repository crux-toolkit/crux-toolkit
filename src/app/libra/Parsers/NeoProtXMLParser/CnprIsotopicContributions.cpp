#include "CnprIsotopicContributions.h"

using namespace std;

void CnprIsotopicContributions::write(FILE* f, int tabs){
  NPRprintTabs(f, tabs);
  fprintf(f, "<isotopic_contributions>\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<contributing_channel.size(); j++) contributing_channel[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</isotopic_contributions>\n");
}