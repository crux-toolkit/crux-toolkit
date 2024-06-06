#include "CnprErrorPoint.h"

using namespace std;

CnprErrorPoint::CnprErrorPoint(){
  error=-1;
  min_prob=-1;
  num_corr=-1;
  num_incorr=-1;
}

void CnprErrorPoint::write(FILE* f, int tabs){
  //required
  string el = "error_point";
  if (error<0) NPRerrMsg(el, "error");
  if (min_prob<0) NPRerrMsg(el, "min_prob");
  if (num_corr<0) NPRerrMsg(el, "num_corr");
  if (num_incorr<0) NPRerrMsg(el, "num_incorr");

  NPRprintTabs(f, tabs);
  fprintf(f, "<error_point");
  fprintf(f, " error=\"%.4lf\"", error);
  fprintf(f, " min_prob=\"%.4lf\"", min_prob);
  fprintf(f, " num_corr=\"%d\"", num_corr);
  fprintf(f, " num_incorr=\"%d\"", num_incorr);
  fprintf(f, "/>\n");

}
