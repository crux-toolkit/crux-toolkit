#include "CnprIndistinguishablePeptide.h"

using namespace std;

CnprIndistinguishablePeptide::CnprIndistinguishablePeptide(){
  charge=-1;
  calc_neutral_pep_mass=-1;
}

void CnprIndistinguishablePeptide::write(FILE* f, int tabs){
  //required
  string el = "indistinguishable_peptide";
  if (peptide_sequence.empty()) NPRerrMsg(el, "peptide_sequence");
  if (charge<0) NPRerrMsg(el, "charge");

  NPRprintTabs(f, tabs);
  fprintf(f, "<indistinguishable_peptide");
  fprintf(f, " peptide_sequence=\"%s\"", peptide_sequence.c_str());
  fprintf(f, " charge=\"%d\"", charge);
  if (calc_neutral_pep_mass>-1)fprintf(f, " calc_neutral_pep_mass=\"%.2lf\"", calc_neutral_pep_mass);
  if (modification_info.empty()) {
    fprintf(f, "/>\n");
    return;
  } else fprintf(f, ">\n");

  int t = tabs;
  if (t>-1) t++;

  for (size_t i = 0; i<modification_info.size(); i++) modification_info[i].write(f, t);

  NPRprintTabs(f, tabs);
  fprintf(f, "</indistinguishable_peptide>\n");

}
