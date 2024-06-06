#include "CnpxModificationInfo.h"

using namespace std;

CnpxModificationInfo::CnpxModificationInfo() {
  modified_peptide.clear();
  mod_cterm_mass = 0;
  mod_nterm_mass = 0;
  mod_aminoacid_mass.clear();
}

CnpxModAminoAcidMass* CnpxModificationInfo::addModAminoAcidMass(int position, double mass){
  CnpxModAminoAcidMass m;
  m.position = position;
  m.mass = mass;
  mod_aminoacid_mass.push_back(m);
  return &mod_aminoacid_mass.back();
}

void CnpxModificationInfo::write(FILE* f, int tabs) {
  size_t i;

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<modification_info");
  if(!modified_peptide.empty()) fprintf(f," modified_peptide=\"%s\"", modified_peptide.c_str());
  if (mod_cterm_mass != 0) fprintf(f, " mod_cterm_mass=\"%.6lf\"", mod_cterm_mass);
  if (mod_nterm_mass != 0) fprintf(f, " mod_nterm_mass=\"%.6lf\"", mod_nterm_mass);
  fprintf(f, ">\n");

  for (i = 0; i<mod_aminoacid_mass.size(); i++) mod_aminoacid_mass[i].write(f,t);
  for (i = 0; i<aminoacid_substitution.size(); i++) aminoacid_substitution[i].write(f, t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</modification_info>\n");
}