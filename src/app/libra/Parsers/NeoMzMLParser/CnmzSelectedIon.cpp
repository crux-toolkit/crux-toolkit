#include "CnmzSelectedIon.h"

using namespace std;

void CnmzSelectedIon::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "selectedIon";

  NMZprintTabs(f, tabs);
  fprintf(f, "<selectedIon>\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroupRef.size(); a++) referenceableParamGroupRef[a].write(f, t, iterative);
  for (size_t a = 0; a<cvParam.size(); a++) cvParam[a].write(f, t, iterative);
  for (size_t a = 0; a<userParam.size(); a++) userParam[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</selectedIon>\n");

}

