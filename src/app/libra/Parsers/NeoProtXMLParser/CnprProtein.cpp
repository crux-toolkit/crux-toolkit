#include "CnprProtein.h"

using namespace std;

CnprProtein::CnprProtein(){
  confidence=-1;
  n_indistinguishable_proteins=-1;
  percent_coverage=-1;
  probability=-1;
  total_number_distinct_peptides=-1;
  total_number_peptides=-1;
}

void CnprProtein::write(FILE* f, int tabs){
  //required
  string el = "protein";
  if (protein_name.empty()) NPRerrMsg(el, "protein_name");
  if (probability<0)NPRerrMsg(el, "probability");
  if (n_indistinguishable_proteins<0)NPRerrMsg(el, "n_indistinguishable_proteins");
  if (group_sibling_id.empty())NPRerrMsg(el, "group_sibling_id");
  if (peptide.empty())NPRerrMsg(el, "peptide");

  int t = tabs;
  if (t>-1) t++;

  NPRprintTabs(f, tabs);
  fprintf(f, "<protein");
  fprintf(f, " protein_name=\"%s\"", protein_name.c_str());
  fprintf(f, " n_indistinguishable_proteins=\"%d\"", n_indistinguishable_proteins);
  fprintf(f, " probability=\"%.4lf\"", probability);
  if (percent_coverage>-1) fprintf(f, " percent_coverage=\"%.1lf\"", percent_coverage);
  if (!unique_stripped_peptides.empty()) fprintf(f, " unique_stripped_peptides=\"%s\"", unique_stripped_peptides.c_str());
  fprintf(f, " group_sibling_id=\"%s\"", group_sibling_id.c_str());
  if (total_number_peptides>-1) fprintf(f, " total_number_peptides=\"%d\"", total_number_peptides);
  if (total_number_distinct_peptides>-1) fprintf(f, " total_number_distinct_peptides=\"%d\"", total_number_distinct_peptides);
  if (!subsuming_protein_entry.empty()) fprintf(f, " subsuming_protein_entry=\"%s\"", subsuming_protein_entry.c_str());
  if (!pct_spectrum_ids.empty()) fprintf(f, " pct_spectrum_ids=\"%s\"", pct_spectrum_ids.c_str());
  if (confidence>-1) fprintf(f, " confidence=\"%.3lf\"", confidence);
  fprintf(f, ">\n");

  size_t j;
  for (j = 0; j<analysis_result.size(); j++) analysis_result[j].write(f, t);
  for(j=0;j<parameter.size();j++) parameter[j].write(f,t);
  for (j = 0; j<annotation.size(); j++) annotation[j].write(f, t);
  for (j = 0; j<indistinguishable_protein.size(); j++) indistinguishable_protein[j].write(f, t);
  for (j = 0; j<peptide.size(); j++) peptide[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</protein>\n");

}



