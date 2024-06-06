#include "CnpxLability.h"

using namespace std;

CnpxLability::CnpxLability(){
  numlosses=0;
  pval=0;
  probability = 0;
  oscore = 0;
  mscore = 0;
  cterm_score = 0;
  nterm_score = 0;
}


void CnpxLability::write(FILE* f) {

  fprintf(f, "<lability numlosses=\"%d\"", numlosses);
  fprintf(f, " pval=\"%.3e\"", pval);
  fprintf(f, " probability=\"%.3e\"", probability);
  fprintf(f, " oscore=\"%.3lf\"", oscore);
  fprintf(f, " mscore=\"%.3lf\"", mscore);
  fprintf(f, " cterm_score=\"%.3lf\"", cterm_score);
  fprintf(f, " nterm_score=\"%.3lf\"", nterm_score);
  fprintf(f, "/>\n");

}