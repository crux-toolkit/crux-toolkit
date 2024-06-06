#include "CnpxPosModelDistribution.h"

using namespace std;

CnpxPosModelDistribution::CnpxPosModelDistribution(){
  type.clear();
  parameter.clear();
}

void CnpxPosModelDistribution::write(FILE* f) {
  size_t i;

  fprintf(f, "<posmodel_distribution");
  if(type.size()>0) fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f,">\n");
  for (i = 0; i < parameter.size(); i++) parameter[i].write(f);
  fprintf(f, "</posmodel_distribution>\n");
}
