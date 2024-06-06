#include "CnmzReferenceableParamGroupRef.h"

using namespace std;

void CnmzReferenceableParamGroupRef::write(FILE* f, int tabs, bool iterative){
  //required
  string el = "referenceableParamGroupRef";

  NMZprintTabs(f, tabs);
  fprintf(f, "<referenceableParamGroupRef");
  fprintf(f, " ref=\"%s\"",ref.c_str());
  fprintf(f,"/>\n");

}

