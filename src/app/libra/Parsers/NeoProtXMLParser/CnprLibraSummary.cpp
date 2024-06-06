#include "CnprLibraSummary.h"

using namespace std;

CnprLibraSummary::CnprLibraSummary(){
  mass_tolerance=0;
  centroiding_preference=-1;
  normalization=-1;
  output_type=-1;
  min_pep_prob=0;
  min_pep_wt=0;
  min_prot_prob=0;
}

void CnprLibraSummary::write(FILE* f, int tabs){
  //required
  string el = "libra_summary";
  if (centroiding_preference<0) NPRerrMsg(el, "centroiding_preference");
  if (normalization<0) NPRerrMsg(el, "normalization");
  if (output_type<0) NPRerrMsg(el, "output_type");

  NPRprintTabs(f, tabs);
  fprintf(f, "<libra_summary");
  if(!version.empty()) fprintf(f, " version=\"%s\"", version.c_str());
  fprintf(f, " mass_tolerance=\"%.3lf\"", mass_tolerance);
  fprintf(f, " centroiding_preference=\"%d\"", centroiding_preference);
  fprintf(f, " normalization=\"%d\"", normalization);
  fprintf(f, " output_type=\"%d\"", output_type);
  if (!channel_code.empty()) fprintf(f, " channel_code=\"%s\"", channel_code.c_str());
  fprintf(f, " min_pep_prob=\"%.4lf\"", min_pep_prob);
  fprintf(f, " min_pep_wt=\"%.4lf\"", min_pep_wt);
  fprintf(f, " min_prot_prob=\"%.4lf\"", min_prot_prob);
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<fragment_masses.size(); j++) fragment_masses[j].write(f, t);
  for (j = 0; j<isotopic_contributions.size(); j++) isotopic_contributions[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</libra_summary>\n");
}