#include "CnpxInterprophetSummary.h"

using namespace std;

CnpxInterprophetSummary::CnpxInterprophetSummary() {
  est_tot_num_correct_pep = 0;
  est_tot_num_correct_psm = 0;
  version.clear();
  options.clear();
}

void CnpxInterprophetSummary::write(FILE* f){
  size_t i;

  fprintf(f, "<interprophet_summary version=\"%s\"", version.c_str());
  fprintf(f, " options=\"%s\"", options.c_str());
  fprintf(f, " est_tot_num_correct_psm=\"%.4lf\"", est_tot_num_correct_psm);
  fprintf(f, " est_tot_num_correct_pep=\"%.4lf\"", est_tot_num_correct_pep);
  fprintf(f, ">\n");

  for (i = 0; i<inputfile.size(); i++) inputfile[i].write(f);
  for (i=0; i<roc_error_data.size();i++) roc_error_data[i].write(f);
  for (i = 0; i<mixturemodel.size(); i++) mixturemodel[i].write(f);
  for (i = 0; i < mixturemodel_distribution.size(); i++) mixturemodel_distribution[i].write(f);

  fprintf(f, "</interprophet_summary>\n");
}