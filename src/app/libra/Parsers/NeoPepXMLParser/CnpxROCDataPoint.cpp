#include "CnpxROCDataPoint.h"

using namespace std;

CnpxROCDataPoint::CnpxROCDataPoint(){
  error = 0;
  min_prob = 0;
  num_corr = 0;
  num_incorr = 0;
  sensitivity = 0;
}

void CnpxROCDataPoint::write(FILE* f){
  fprintf(f, "<roc_data_point min_prob=\"%.4lf\"", min_prob);
  fprintf(f, " sensitivity=\"%.4lf\"", sensitivity);
  fprintf(f, " error=\"%.4lf\"", error);
  fprintf(f, " num_corr=\"%u\"", num_corr);
  fprintf(f, " num_incorr=\"%u\"", num_incorr);
  fprintf(f, "/>\n");
}