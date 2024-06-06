#include "CnpxSearchScoreSummary.h"

using namespace std;

CnpxSearchScoreSummary::CnpxSearchScoreSummary() {
  parameter.clear();
  active = false;
}

CnpxSearchScoreSummary::CnpxSearchScoreSummary(bool b) {
  parameter.clear();
  active = b;
}

bool CnpxSearchScoreSummary::present() {
  return active;
}

void CnpxSearchScoreSummary::write(FILE* f) {
  size_t i;

  fprintf(f, "<search_score_summary>\n");
  for (i = 0; i < parameter.size(); i++) parameter[i].write(f);
  fprintf(f, "</search_score_summary>\n");
}