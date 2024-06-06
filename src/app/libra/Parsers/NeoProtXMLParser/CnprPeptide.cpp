#include "CnprPeptide.h"

using namespace std;

CnprPeptide::CnprPeptide(){
  charge=-1;
  initial_probability=-1;
  nsp_adjusted_probability=-1;
  fpkm_adjusted_probability=-1;
  ni_adjusted_probability=-1;
  exp_sibling_ion_instances=-1;
  exp_sibling_ion_bin=-1;
  exp_tot_instances=-1;
  weight=1;
  n_enzymatic_termini=-1;
  n_sibling_peptides=-1;
  n_sibling_peptides_bin=0;
  max_fpkm=-1;
  fpkm_bin=0;
  n_instances=-1;
  calc_neutral_pep_mass=-1;
}

void CnprPeptide::write(FILE* f, int tabs){
  //required
  string el = "peptide";
  if (peptide_sequence.empty()) NPRerrMsg(el, "peptide_sequence");
  if (charge<0)NPRerrMsg(el, "charge");
  if (initial_probability<0)NPRerrMsg(el, "initial_probability");
  if (is_nondegenerate_evidence.empty())NPRerrMsg(el, "is_nondegenerate_evidence");
  if (n_enzymatic_termini<0)NPRerrMsg(el, "n_enzymatic_termini");
  if (n_instances<0)NPRerrMsg(el, "n_instances");
  if (is_contributing_evidence.empty())NPRerrMsg(el, "is_contributing_evidence");

  NPRprintTabs(f, tabs);
  fprintf(f, "<peptide");
  fprintf(f, " peptide_sequence=\"%s\"", peptide_sequence.c_str());
  fprintf(f, " charge=\"%d\"", charge);
  fprintf(f, " initial_probability=\"%.4lf\"", initial_probability);
  if (nsp_adjusted_probability>-1) fprintf(f, " nsp_adjusted_probability=\"%.4lf\"", nsp_adjusted_probability);
  if (fpkm_adjusted_probability>-1) fprintf(f, " fpkm_adjusted_probability=\"%.4lf\"", fpkm_adjusted_probability);
  if (ni_adjusted_probability>-1) fprintf(f, " ni_adjusted_probability=\"%.4lf\"", ni_adjusted_probability);
  if (exp_sibling_ion_instances>-1) fprintf(f, " exp_sibling_ion_instances=\"%.4lf\"", exp_sibling_ion_instances);
  if (exp_sibling_ion_bin>-1) fprintf(f, " exp_sibling_ion_bin=\"%.4lf\"", exp_sibling_ion_bin);
  if (exp_tot_instances>-1) fprintf(f, " exp_tot_instances=\"%.4lf\"", exp_tot_instances);
  if (!peptide_group_designator.empty()) fprintf(f, " peptide_group_designator=\"%s\"", peptide_group_designator.c_str());
  fprintf(f, " weight=\"%.2lf\"", weight);
  fprintf(f, " is_nondegenerate_evidence=\"%s\"", is_nondegenerate_evidence.c_str());
  fprintf(f, " n_enzymatic_termini=\"%d\"", n_enzymatic_termini);
  if (n_sibling_peptides>-1) fprintf(f, " n_sibling_peptides=\"%.2lf\"", n_sibling_peptides);
  if (n_sibling_peptides_bin>-1) fprintf(f, " n_sibling_peptides_bin=\"%d\"", n_sibling_peptides_bin);
  if (max_fpkm>-1) fprintf(f, " max_fpkm=\"%.2lf\"", max_fpkm);
  if (fpkm_bin>0) fprintf(f, " fpkm_bin=\"%d\"", fpkm_bin);
  fprintf(f, " n_instances=\"%d\"", n_instances);
  if (calc_neutral_pep_mass>-1) fprintf(f, " calc_neutral_pep_mass=\"%.2lf\"", calc_neutral_pep_mass);
  fprintf(f, " is_contributing_evidence=\"%s\"", is_contributing_evidence.c_str());

  if (parameter.empty() && modification_info.empty() && peptide_parent_protein.empty() && indistinguishable_peptide.empty()){
    fprintf(f, "/>\n");
    return;
  } else fprintf(f, ">\n");
  
  int t = tabs;
  if (t>-1) t++;

  size_t j;
  for (j = 0; j<parameter.size(); j++) parameter[j].write(f, t);
  for (j = 0; j<modification_info.size(); j++) modification_info[j].write(f, t);
  for (j = 0; j<peptide_parent_protein.size(); j++) peptide_parent_protein[j].write(f, t);
  for (j = 0; j<indistinguishable_peptide.size(); j++) indistinguishable_peptide[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</peptide>\n");
}

