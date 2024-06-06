#include "CnpxLibraResult.h"

using namespace std;

CnpxLibraResult::CnpxLibraResult() {
  active = false;
}

CnpxLibraResult::CnpxLibraResult(bool b) {
  active = b;
}

bool CnpxLibraResult::present() {
  return active;
}

void CnpxLibraResult::write(FILE* f, int tabs) {

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<libra_result>\n");

  for(size_t i=0;i<intensity.size();i++) intensity[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</libra_result>\n");
}