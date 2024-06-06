#include "CnpxPTMProphetResult.h"

using namespace std;

CnpxPTMProphetResult::CnpxPTMProphetResult() {
  active = false;
}

CnpxPTMProphetResult::CnpxPTMProphetResult(bool b) {
  active = b;
}

bool CnpxPTMProphetResult::present() {
  return active;
}

void CnpxPTMProphetResult::write(FILE* f) {
  fprintf(f, "<ptmprophet_result ptm=\"%s\"", ptm.c_str());
  if (!prior.empty()) fprintf(f, " prior=\"%s\"", prior.c_str());
  if (!ptm_peptide.empty()) fprintf(f, " ptm_peptide=\"%s\"", ptm_peptide.c_str());
  fprintf(f, ">\n");

  size_t i;
  for(i=0;i<parameter.size();i++) parameter[i].write(f);
  for (i = 0; i<lability.size(); i++) lability[i].write(f);
  for (i = 0; i<mod_terminal_probability.size(); i++) mod_terminal_probability[i].write(f);
  for (i = 0; i<mod_amino_acid_probability.size(); i++) mod_amino_acid_probability[i].write(f);

  fprintf(f, "</ptmprophet_result>\n");
}