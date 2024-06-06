#include "CnprProteinProphetDetails.h"

using namespace std;

void CnprProteinProphetDetails::write(FILE* f, int tabs){
  //required
  string el = "proteinprophet_details";
  if (occam_flag.empty()) NPRerrMsg(el, "occam_flag");
  if (groups_flag.empty()) NPRerrMsg(el, "groups_flag");
  if (degen_flag.empty()) NPRerrMsg(el, "degen_flag");
  if (nsp_flag.empty()) NPRerrMsg(el, "nsp_flag");
  if (fpkm_flag.empty()) NPRerrMsg(el, "fpkm_flag");
  if (initial_peptide_wt_iters.empty()) NPRerrMsg(el, "initial_peptide_wt_iters");
  if (nsp_distribution_iters.empty()) NPRerrMsg(el, "nsp_distribution_iters");
  if (final_peptide_wt_iters.empty()) NPRerrMsg(el, "final_peptide_wt_iters");
  //if (protein_summary_data_filter.empty())  NPRerrMsg(el, "protein_summary_data_filter");

  NPRprintTabs(f, tabs);
  fprintf(f, "<proteinprophet_details");
  fprintf(f, " occam_flag=\"%s\"", occam_flag.c_str());
  fprintf(f, " groups_flag=\"%s\"", groups_flag.c_str());
  fprintf(f, " degen_flag=\"%s\"", degen_flag.c_str());
  fprintf(f, " nsp_flag=\"%s\"", nsp_flag.c_str());
  fprintf(f, " fpkm_flag=\"%s\"", fpkm_flag.c_str());
  fprintf(f, " initial_peptide_wt_iters=\"%s\"", initial_peptide_wt_iters.c_str());
  fprintf(f, " nsp_distribution_iters=\"%s\"", nsp_distribution_iters.c_str());
  fprintf(f, " final_peptide_wt_iters=\"%s\"", final_peptide_wt_iters.c_str());
  if (!run_options.empty()) fprintf(f, " run_options=\"%s\"", run_options.c_str());
  fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  nsp_information.write(f,t);

  size_t j;
  //for (j = 0; j<fpkm_information.size(); j++) fpkm_information[j].write(f, t);
  //for (j = 0; j<ni_information.size(); j++) ni_information[j].write(f, t);
  for (j = 0; j<protein_summary_data_filter.size(); j++) protein_summary_data_filter[j].write(f, t);
  for (j = 0; j<error_point.size(); j++) error_point[j].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</proteinprophet_details>\n");

}
