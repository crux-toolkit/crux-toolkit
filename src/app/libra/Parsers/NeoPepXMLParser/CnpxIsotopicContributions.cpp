#include "CnpxIsotopicContributions.h"

using namespace std;

void CnpxIsotopicContributions::write(FILE* f, int tabs){
  size_t i;

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<isotopic_contributions>\n");

  for (i = 0; i<contributing_channel.size(); i++) contributing_channel[i].write(f, t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</isotopic_contributions>\n");
}