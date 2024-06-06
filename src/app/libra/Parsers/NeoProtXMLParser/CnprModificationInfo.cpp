#include "CnprModificationInfo.h"

using namespace std;

CnprModificationInfo::CnprModificationInfo(){
  mod_cterm_mass=0;
  mod_nterm_mass=0;
}

void CnprModificationInfo::write(FILE* f, int tabs){
  //nothing required

  NPRprintTabs(f, tabs);
  fprintf(f, "<modification_info");
  if (mod_nterm_mass!=0)fprintf(f, " mod_nterm_mass=\"%.6lf\"", mod_nterm_mass);
  if (mod_cterm_mass!=0)fprintf(f, " mod_cterm_mass=\"%.6lf\"", mod_cterm_mass);
  if (!modified_peptide.empty())fprintf(f, " modified_peptide=\"%s\"", modified_peptide.c_str());
  if(mod_aminoacid_mass.empty()) {
    fprintf(f, "/>\n");
    return;
  } else fprintf(f,">\n");

  int t = tabs;
  if (t>-1) t++;

  for(size_t i=0;i<mod_aminoacid_mass.size();i++) mod_aminoacid_mass[i].write(f,t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</modification_info>");

}
