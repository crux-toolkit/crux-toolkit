#include "CnpxLibraSummary.h"

using namespace std;

CnpxLibraSummary::CnpxLibraSummary(){
  mass_tolerance=0;
  centroiding_preference=0;
  normalization=0;
  output_type=0;
}

void CnpxLibraSummary::write(FILE* f, int tabs){
  size_t i;

  string el = "libra_summary";

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<libra_summary");
  fprintf(f, " mass_tolerance=\"%3lf\"", mass_tolerance);
  fprintf(f, " centroiding_preference=\"%d\"", centroiding_preference);
  fprintf(f, " normalization=\"%d\"", normalization);
  fprintf(f, " output_type=\"%d\"", output_type);
  fprintf(f, " channel_code=\"%s\"", channel_code.c_str());
  fprintf(f, ">\n");

  for (i = 0; i<fragment_masses.size(); i++) fragment_masses[i].write(f,t);
  for (i = 0; i<isotopic_contributions.size(); i++) isotopic_contributions[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</libra_summary>\n");
}