#include "CnpxAminoAcidSubstitution.h"

using namespace std;

CnpxAminoAcidSubstitution::CnpxAminoAcidSubstitution() {
  position = 0;
  num_tol_term = 0;
}

void CnpxAminoAcidSubstitution::write(FILE* f, int tabs) {
  string el = "aminoacid_substitution";

  NPXprintTabs(f, tabs);
  fprintf(f, "<aminoacid_substitution");
  if(position!=0) fprintf(f, " position=\"%d\"", position);
  if (!orig_aa.empty()) fprintf(f, " orig_aa=\"%s\"", orig_aa.c_str());
  fprintf(f, " num_tol_term=\"%d\"", num_tol_term);
  if (!peptide_prev_aa.empty()) fprintf(f, " peptide_prev_aa=\"%s\"", peptide_prev_aa.c_str());
  if (!peptide_next_aa.empty()) fprintf(f, " peptide_next_aa=\"%s\"", peptide_next_aa.c_str());
  fprintf(f, "/>\n");
}
