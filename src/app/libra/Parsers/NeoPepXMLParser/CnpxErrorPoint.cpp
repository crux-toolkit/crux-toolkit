#include "CnpxErrorPoint.h"

using namespace std;

CnpxErrorPoint::CnpxErrorPoint(){
  error = 0;
  min_prob = 0;
  num_corr = 0;
  num_incorr = 0;
}

void CnpxErrorPoint::write(FILE* f){
  fprintf(f, "<error_point error=\"%.4lf\"", error);
  fprintf(f, " min_prob=\"%.4lf\"", min_prob);
  fprintf(f, " num_corr=\"%u\"", num_corr);
  fprintf(f, " num_incorr=\"%u\"", num_incorr);
  fprintf(f, "/>\n");
}