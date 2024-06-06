#include "CnprLibraResult.h"

using namespace std;

CnprLibraResult::CnprLibraResult(){
  number = -1;
}

void CnprLibraResult::write(FILE* f, int tabs){
  //required
  string el = "libra_result";
  if (number<0) NPRerrMsg(el, "number");

  NPRprintTabs(f, tabs);
  fprintf(f, "<libra_result");
  fprintf(f, " number=\"%d\"", number);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<intensity.size(); j++) intensity[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</libra_result>\n");

}
