#include "CnpxPeptideprophetSummary.h"

using namespace std;

CnpxPeptideprophetSummary::CnpxPeptideprophetSummary() {
  double est_tot_num_correct=0;
  double min_prob=0;
  author.clear();
  options.clear();
  type.clear();
  version.clear();
}

void CnpxPeptideprophetSummary::write(FILE* f) {
  size_t i;

  fprintf(f, "<peptideprophet_summary version=\"%s\"", version.c_str());
  fprintf(f, " author=\"%s\"", author.c_str());
  fprintf(f, " type=\"%s\"", type.c_str());
  fprintf(f, " min_prob=\"%.2lf\"", min_prob);
  fprintf(f, " options=\"%s\"", options.c_str());
  fprintf(f, " est_tot_num_correct=\"%.4lf\"", est_tot_num_correct);
  fprintf(f, ">\n");

  for (i = 0; i < inputfile.size(); i++) inputfile[i].write(f);
  for (i = 0; i < roc_error_data.size(); i++) roc_error_data[i].write(f);
  for (i = 0; i < distribution_point.size(); i++) distribution_point[i].write(f);
  for (i = 0; i < mixture_model.size(); i++) mixture_model[i].write(f);

  fprintf(f, "</peptideprophet_summary>\n");
}