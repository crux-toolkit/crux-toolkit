#include "CnpxMixtureModelDistribution.h"

using namespace std;

CnpxMixtureModelDistribution::CnpxMixtureModelDistribution() {
  name.clear();
  negmodel_distribution.clear();
  posmodel_distribution.clear();
}

void CnpxMixtureModelDistribution::write(FILE* f) {
  size_t i;

  fprintf(f, "<mixturemodel_distribution name=\"%s\">\n", name.c_str());
  for (i = 0; i < posmodel_distribution.size(); i++) posmodel_distribution[i].write(f);
  for (i = 0; i < negmodel_distribution.size(); i++) negmodel_distribution[i].write(f);
  fprintf(f, "</mixturemodel_distribution>\n");
}
