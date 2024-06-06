#include "CnprProteinSummary.h"

using namespace std;

void CnprProteinSummary::write(FILE* f, int tabs){
  NPRprintTabs(f,tabs);
  fprintf(f,"<protein_summary xmlns=\"%s\"", xmlns.c_str());
  fprintf(f," xmlns:xsi=\"%s\"", xmlns_xsi.c_str());
  fprintf(f," xsi:schemaLocation=\"%s\"", xsi_schemaLocation.c_str());
  fprintf(f," summary_xml=\"%s\">\n", summary_xml.c_str());

  int t=tabs;
  if(t>-1) t++;

  size_t j;
  protein_summary_header.write(f,t);
  for (j = 0; j<analysis_summary.size(); j++) analysis_summary[j].write(f, t);
  dataset_derivation.write(f,t);
  for(j=0;j<protein_group.size();j++) protein_group[j].write(f,t);

  NPRprintTabs(f,tabs);
  fprintf(f,"</protein_summary>\n");
}