#include "CnpxModAminoAcidMass.h"

using namespace std;

CnpxModAminoAcidMass::CnpxModAminoAcidMass() {
  id.clear();
  mass=0;
  position=0;
  source.clear();
  staticMass=0;
  variable=0;
}

void CnpxModAminoAcidMass::write(FILE* f, int tabs) {
  string el = "mod_aminoacid_mass";
  if (position==0) NPXerrMsg(el, "position");
  if (mass == 0) NPXerrMsg(el, "mass");

  NPXprintTabs(f, tabs);
  fprintf(f, "<mod_aminoacid_mass position=\"%d\"", position);
  fprintf(f, " mass=\"%.6lf\"", mass);
  if (staticMass != 0) fprintf(f, " static=\"%.6lf\"", staticMass);
  if (variable != 0) fprintf(f, " variable=\"%.6lf\"", variable);
  if (source.size() > 0) fprintf(f, " source=\"%s\"", source.c_str());
  if (id.size() > 0) fprintf(f, " id=\"%s\"", id.c_str());
  fprintf(f, "/>\n");
}
