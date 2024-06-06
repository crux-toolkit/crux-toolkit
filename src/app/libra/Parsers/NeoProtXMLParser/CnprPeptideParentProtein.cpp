#include "CnprPeptideParentProtein.h"

using namespace std;

void CnprPeptideParentProtein::write(FILE* f, int tabs){
  //required
  string el = "peptide_parent_protein";
  if (protein_name.empty()) NPRerrMsg(el, "protein_name");

  NPRprintTabs(f, tabs);
  fprintf(f, "<peptide_parent_protein");
  fprintf(f, " protein_name=\"%s\"", protein_name.c_str());
  fprintf(f, "/>\n");

}
