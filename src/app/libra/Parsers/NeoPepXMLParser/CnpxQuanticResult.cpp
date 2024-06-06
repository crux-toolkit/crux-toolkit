#include "CnpxQuanticResult.h"

using namespace std;

CnpxQuanticResult::CnpxQuanticResult() {
  antic=-1;
  active = false;
}

CnpxQuanticResult::CnpxQuanticResult(bool b) {
  antic=-1;
  active = b;
}

bool CnpxQuanticResult::present() {
  return active;
}

void CnpxQuanticResult::write(FILE* f) {
  fprintf(f, "<quantic_result antic=\"%.2lf\"/>\n", antic);
}