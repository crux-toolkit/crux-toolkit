#include "CnpxEnzymaticSearchConstraint.h"

using namespace std;

CnpxEnzymaticSearchConstraint::CnpxEnzymaticSearchConstraint(){
  enzyme.clear();
  max_num_internal_cleavages=0;
  min_number_termini=0;
}

void CnpxEnzymaticSearchConstraint::write(FILE* f){

  fprintf(f, "<enzymatic_search_constraint enzyme=\"%s\"", enzyme.c_str());
  fprintf(f, " max_num_internal_cleavages=\"%d\"", max_num_internal_cleavages);
  fprintf(f, " min_number_termini=\"%d\"", min_number_termini);
  fprintf(f, "/>\n");

}
