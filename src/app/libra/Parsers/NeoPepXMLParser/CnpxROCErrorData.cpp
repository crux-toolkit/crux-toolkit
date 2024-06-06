#include "CnpxROCErrorData.h"

using namespace std;

CnpxROCErrorData::CnpxROCErrorData(){
  charge.clear();
  charge_est_correct=0;
}

void CnpxROCErrorData::write(FILE* f){
  size_t i;

  fprintf(f, "<roc_error_data charge=\"%s\"", charge.c_str());
  if(charge_est_correct!=0) fprintf(f, " charge_est_correct=\"%.1lf\"", charge_est_correct);
  fprintf(f, ">\n");

  for (i = 0; i<roc_data_point.size(); i++) roc_data_point[i].write(f);
  for (i = 0; i<error_point.size(); i++) error_point[i].write(f);

  fprintf(f, "</roc_error_data>\n");
}