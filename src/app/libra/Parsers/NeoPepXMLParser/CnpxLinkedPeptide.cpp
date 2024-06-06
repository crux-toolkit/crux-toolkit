#include "CnpxLinkedPeptide.h"

using namespace std;

CnpxLinkedPeptide::CnpxLinkedPeptide(){
  designation.clear();
  peptide.clear();
  peptide_prev_aa.clear();
  peptide_next_aa.clear();
  protein.clear();
  num_tot_proteins = 0;
  calc_neutral_pep_mass = 0;
  protein_link_pos_a=0;
  complement_mass=0;
}

void CnpxLinkedPeptide::write(FILE* f, int tabs){
  size_t i;

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<linked_peptide peptide=\"%s\"", peptide.c_str());
  if (peptide_prev_aa.size()>0) fprintf(f, " peptide_prev_aa=\"%s\"", peptide_prev_aa.c_str());
  if (peptide_next_aa.size()>0) fprintf(f, " peptide_next_aa=\"%s\"", peptide_next_aa.c_str());
  fprintf(f, " protein=\"%s\"", protein.c_str());
  if (peptide_start_pos>0) fprintf(f, " peptide_start_pos=\"%d\"", peptide_start_pos);
  if (protein_link_pos_a>0) fprintf(f, " protein_link_pos_a=\"%d\"", protein_link_pos_a);
  fprintf(f, " num_tot_proteins=\"%d\"", num_tot_proteins);
  fprintf(f, " calc_neutral_pep_mass=\"%.6lf\"", calc_neutral_pep_mass);
  fprintf(f, " complement_mass=\"%.6lf\"", complement_mass);
  if (designation.size()>0) fprintf(f, " designation=\"%s\"", designation.c_str());
  fprintf(f, ">\n");

  for (i = 0; i < alternative_protein.size(); i++) alternative_protein[i].write(f,t);
  for (i = 0; i < modification_info.size(); i++) modification_info[i].write(f, t);
  for (i = 0; i < xlink_score.size(); i++) xlink_score[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</linked_peptide>\n");

}
