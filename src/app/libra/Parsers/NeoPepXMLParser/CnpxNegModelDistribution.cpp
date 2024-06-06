#include "CnpxNegModelDistribution.h"

using namespace std;

CnpxNegModelDistribution::CnpxNegModelDistribution() {
  type.clear();
  parameter.clear();
}

void CnpxNegModelDistribution::write(FILE* f) {
  size_t i;

  fprintf(f, "<negmodel_distribution");
  if (type.size() > 0) fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f, ">\n");
  for (i = 0; i < parameter.size(); i++) parameter[i].write(f);
  fprintf(f, "</negmodel_distribution>\n");
}
