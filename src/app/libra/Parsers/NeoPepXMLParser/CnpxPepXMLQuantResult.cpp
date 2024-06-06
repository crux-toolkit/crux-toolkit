#include "CnpxPepXMLQuantResult.h"

using namespace std;

CnpxPepXMLQuantResult::CnpxPepXMLQuantResult() {
  area = 0;
  retention_time_sec = 0;
  active = false;
}

CnpxPepXMLQuantResult::CnpxPepXMLQuantResult(bool b) {
  area = 0;
  retention_time_sec = 0;
  active = b;
}

bool CnpxPepXMLQuantResult::present() {
  return active;
}

void CnpxPepXMLQuantResult::write(FILE* f) {
  fprintf(f, "<pepxmlquant_result area=\"%.2lf\" retention_time_sec=\"%.4lf\"", area, retention_time_sec);
  fprintf(f, ">\n");
  if (search_score_summary.present()) search_score_summary.write(f);
  fprintf(f, "</pepxmlquant_result>\n");
}