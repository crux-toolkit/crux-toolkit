#include "CnprProteinSummaryHeader.h"

using namespace std;

CnprProteinSummaryHeader::CnprProteinSummaryHeader(){
  initial_min_peptide_pro=-1;
  min_peptide_probability=-1;
  min_peptide_weight=-1;
  num_input_1_spectra=-1;
  num_input_2_spectra=-1;
  num_input_3_spectra=-1;
  num_input_4_spectra=-1;
  num_input_5_spectra=-1;
  num_predicted_correct_prots=-1;
  total_no_spectrum_ids=-1;
}

void CnprProteinSummaryHeader::write(FILE* f, int tabs){
  //required
  string el="protein_summary_header";
  if(initial_min_peptide_pro<0) NPRerrMsg(el,"initial_min_peptide_pro");
  if (min_peptide_probability<0)NPRerrMsg(el, "min_peptide_probability");
  if (min_peptide_weight<0)NPRerrMsg(el, "min_peptide_weight");
  if (num_input_1_spectra<0)NPRerrMsg(el, "num_input_1_spectra");
  if (num_input_2_spectra<0)NPRerrMsg(el, "num_input_2_spectra");
  if (num_input_3_spectra<0)NPRerrMsg(el, "num_input_3_spectra");
  if (num_input_4_spectra<0)NPRerrMsg(el, "num_input_4_spectra");
  if (num_input_5_spectra<0)NPRerrMsg(el, "num_input_5_spectra");
  if (num_predicted_correct_prots<0)NPRerrMsg(el, "num_predicted_correct_prots");
  if (reference_database.empty())NPRerrMsg(el, "reference_database");
  if (residue_substitution_list.empty())NPRerrMsg(el, "residue_substitution_list");
  if (sample_enzyme.empty())NPRerrMsg(el, "sample_enzyme");
  if (source_files.empty())NPRerrMsg(el, "source_files");
  if (source_files_alt.empty())NPRerrMsg(el, "source_files_alt");

  NPRprintTabs(f, tabs);
  fprintf(f, "<protein_summary_header");
  fprintf(f, " reference_database=\"%s\"",reference_database.c_str());
  fprintf(f, " residue_substitution_list=\"%s\"", residue_substitution_list.c_str());
  fprintf(f, " source_files=\"%s\"", source_files.c_str());
  fprintf(f, " source_files_alt=\"%s\"", source_files_alt.c_str());
  fprintf(f, " initial_min_peptide_pro=\"%.2lf\"", initial_min_peptide_pro);
  fprintf(f, " min_peptide_probability=\"%.2lf\"", min_peptide_probability);
  fprintf(f, " min_peptide_weight=\"%.2lf\"", min_peptide_weight);
  fprintf(f, " num_input_1_spectra=\"%d\"", num_input_1_spectra);
  fprintf(f, " num_input_2_spectra=\"%d\"", num_input_2_spectra);
  fprintf(f, " num_input_3_spectra=\"%d\"", num_input_3_spectra);
  fprintf(f, " num_input_4_spectra=\"%d\"", num_input_4_spectra);
  fprintf(f, " num_input_5_spectra=\"%d\"", num_input_5_spectra);
  fprintf(f, " num_predicted_correct_prots=\"%.1lf\"", num_predicted_correct_prots);
  fprintf(f, " sample_enzyme=\"%s\"",sample_enzyme.c_str());

  if (!organism.empty()) fprintf(f, " organism=\"%s\"", organism.c_str());
  if (!source_file_xtn.empty()) fprintf(f, " source_file_xtn=\"%s\"", source_file_xtn.c_str());
  if (total_no_spectrum_ids>-1) fprintf(f, " total_no_spectrum_ids=\"%.1lf\"", total_no_spectrum_ids);
  if (!win_cyg_reference_database.empty()) fprintf(f, " win_cyg_reference_database=\"%s\"", win_cyg_reference_database.c_str());

  fprintf(f,">\n");

  int t = tabs;
  if (t>-1) t++;

  program_details.write(f,t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</protein_summary_header>\n");
}
