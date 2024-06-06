#include "CnpxMixture_Model.h"

using namespace std;

CnpxMixture_Model::CnpxMixture_Model(){
  comments.clear();
  est_tot_correct=0;
  num_iterations=0;
  precursor_ion_charge=0;
  prior_probability=0;
  tot_num_spectra=0;

  mixturemodel_distribution.clear();
  mixturemodel.clear();
}

void CnpxMixture_Model::write(FILE* f) {
  size_t i;

  fprintf(f, "<mixture_model precursor_ion_charge=\"%d\"", precursor_ion_charge);
  fprintf(f, " comments=\"%s\"",comments.c_str());
  fprintf(f, " prior_probability=\"%.3lf\"",prior_probability);
  fprintf(f, " est_tot_correct=\"%.1lf\"",est_tot_correct);
  fprintf(f, " tot_num_spectra=\"%d\"",tot_num_spectra);
  fprintf(f, " num_iterations=\"%d\"",num_iterations);
  fprintf(f, ">\n");
  for (i = 0; i < mixturemodel_distribution.size(); i++) mixturemodel_distribution[i].write(f);
  for (i = 0; i < mixturemodel.size(); i++) mixturemodel[i].write(f);
  fprintf(f, "</mixture_model>\n");
}