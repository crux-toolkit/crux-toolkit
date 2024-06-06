#include "CnpxSampleEnzyme.h"

using namespace std;

CnpxSampleEnzyme::CnpxSampleEnzyme(){
  description.clear();
  fidelity.clear();
  independent=true;
  name.clear();

  specificity.clear();
}

void CnpxSampleEnzyme::write(FILE* f, int tabs){
  size_t i;

  string el = "sample_enzyme";
  if (name.empty()) NPXerrMsg(el, "name");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);
  
  fprintf(f, "<sample_enzyme name=\"%s\"", name.c_str());
  if(!description.empty()) fprintf(f, " description=\"%s\"", description.c_str());
  if (!fidelity.empty()) fprintf(f, " fidelity=\"%s\"", fidelity.c_str());
  if (!independent) fprintf(f, " independent=\"0");
  fprintf(f, ">\n");

  for(i=0;i<specificity.size();i++) specificity[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f,"</sample_enzyme>\n");

}