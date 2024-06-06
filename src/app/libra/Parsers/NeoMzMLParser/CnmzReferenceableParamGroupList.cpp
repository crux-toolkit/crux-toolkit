#include "CnmzReferenceableParamGroupList.h"

using namespace std;

void CnmzReferenceableParamGroupList::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "referenceableParamGroupList";
  if (referenceableParamGroup.empty()) NMZerrMsg(el, "referenceableParamGroup");
  if (count != referenceableParamGroup.size()) NMZerrMsg(el, "count match number of referenceableParamGroup");

  NMZprintTabs(f, tabs);
  fprintf(f, "<referenceableParamGroupList");
  fprintf(f, " count=\"%d\"",count);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;
  for (size_t a = 0; a<referenceableParamGroup.size(); a++) referenceableParamGroup[a].write(f, t, iterative);

  NMZprintTabs(f, tabs);
  fprintf(f, "</referenceableParamGroupList>\n");

}

