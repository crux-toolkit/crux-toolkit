#include "CnpxPTMProphetSummary.h"

using namespace std;

CnpxPTMProphetSummary::CnpxPTMProphetSummary() {
  frag_ppm_tol.clear();
  min_o.clear();
  min_o_factors.clear();
  mod_string.clear();
  options.clear();
  version.clear();
  version.clear();
  options.clear();
}

void CnpxPTMProphetSummary::write(FILE* f){
  size_t i;

  fprintf(f, "<ptmprophet_summary version=\"%s\"", version.c_str());
  fprintf(f, " options=\"%s\"", options.c_str());
  if(!options.empty()) fprintf(f, " options=\"%s\"", options.c_str());
  if (!mod_string.empty()) fprintf(f, " mod_string=\"%s\"", mod_string.c_str());
  if (!min_o.empty()) fprintf(f, " min_o=\"%s\"", min_o.c_str());
  if (!min_o_factors.empty()) fprintf(f, " min_o_factors=\"%s\"", min_o_factors.c_str());
  if (!frag_ppm_tol.empty()) fprintf(f, " frag_ppm_tol=\"%s\"", frag_ppm_tol.c_str());
  fprintf(f, ">\n");

  for (i = 0; i<roc_error_data.size(); i++) roc_error_data[i].write(f);
  for (i = 0; i<inputfile.size(); i++) inputfile[i].write(f);
  for (i = 0; i<mixturemodel.size(); i++) mixturemodel[i].write(f);

  fprintf(f, "</interprophet_summary>\n");
}