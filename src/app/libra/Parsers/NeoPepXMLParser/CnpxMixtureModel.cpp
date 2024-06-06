#include "CnpxMixtureModel.h"

using namespace std;

CnpxMixtureModel::CnpxMixtureModel(){
  name.clear();
  neg_bandwidth = 0;
  pos_bandwidth = 0;
  point.clear();
}

void CnpxMixtureModel::write(FILE* f){
  size_t i;

  fprintf(f, "<mixturemodel name=\"%s\"", name.c_str());
  fprintf(f, " pos_bandwidth=\"%.6f\"", pos_bandwidth);
  fprintf(f, " neg_bandwidth=\"%.6f\"", pos_bandwidth);
  fprintf(f, ">\n");

  for (i = 0; i<point.size(); i++) point[i].write(f);
  for (i = 0; i < bin.size(); i++) bin[i].write(f);

  fprintf(f, "</mixturemodel>\n");
}