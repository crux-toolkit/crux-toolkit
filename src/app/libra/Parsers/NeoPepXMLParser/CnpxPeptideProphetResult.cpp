#include "CnpxPeptideProphetResult.h"

using namespace std;

CnpxPeptideProphetResult::CnpxPeptideProphetResult() {
  all_ntt_prob.clear();
  analysis.clear();
  probability = -1;
  pep1_probability=-1;
  pep2_probability=-1;
  active = false;
}

CnpxPeptideProphetResult::CnpxPeptideProphetResult(bool b) {
  all_ntt_prob.clear();
  analysis.clear();
  probability = -1;
  pep1_probability = -1;
  pep2_probability = -1;
  active = b;
}

bool CnpxPeptideProphetResult::present() {
  return active;
}

void CnpxPeptideProphetResult::write(FILE* f) {
  fprintf(f, "<peptideprophet_result probability=\"%.6lf\"", probability);
  if (all_ntt_prob.size() > 0) fprintf(f, " all_ntt_prob=\"%s\"", all_ntt_prob.c_str());
  if (pep1_probability>=0) fprintf(f, " pep1_probability=\"%.6lf\"", pep1_probability);
  if (pep2_probability >= 0) fprintf(f, " pep2_probability=\"%.6lf\"", pep2_probability);
  if (analysis.size() > 0) fprintf(f, " analysis=\"%s\"", analysis.c_str());
  fprintf(f, ">\n");

  if (search_score_summary.present()) search_score_summary.write(f);

  fprintf(f, "</peptideprophet_result>\n");
}