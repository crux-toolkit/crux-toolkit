#include "CnpxSearchScore.h"

using namespace std;

void CnpxSearchScore::write(FILE* f, int tabs) {
  NPXprintTabs(f, tabs);
  fprintf(f, "<search_score name=\"%s\" value=\"%s\"/>\n", name.c_str(), value.c_str());
}