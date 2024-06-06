#include "CnpxInterprophetResult.h"

using namespace std;

CnpxInterprophetResult::CnpxInterprophetResult() {
  all_ntt_prob.clear();
  probability = 0;
  active = false;
}

CnpxInterprophetResult::CnpxInterprophetResult(bool b) {
  all_ntt_prob.clear();
  probability = 0;
  active = b;
}

bool CnpxInterprophetResult::present() {
  return active;
}

void CnpxInterprophetResult::write(FILE* f) {
  fprintf(f, "<interprophet_result probability=\"%.6lf\"", probability);
  if (all_ntt_prob.size() > 0) fprintf(f, " all_ntt_prob=\"%s\"", all_ntt_prob.c_str());
  fprintf(f, ">\n");

  if (search_score_summary.present()) search_score_summary.write(f);

  fprintf(f, "</interprophet_result>\n");
}