#include "CnpxXpressLabelFreeSummary.h"

using namespace std;

CnpxXpressLabelFreeSummary::CnpxXpressLabelFreeSummary() {

}

void CnpxXpressLabelFreeSummary::write(FILE* f){
  fprintf(f, "<xpresslabelfree_summary version=\"%s\"", version.c_str());
  fprintf(f, " author=\"%s\"", author.c_str());
  fprintf(f, " masstol=\"%s\"", masstol.c_str());
  fprintf(f, " ppmtol=\"%s\"", ppmtol.c_str());
  fprintf(f, " min_num_chromatogram_points=\"%s\"", min_num_chromatogram_points.c_str());
  fprintf(f, " min_num_isotope_peaks=\"%s\"", min_num_isotope_peaks.c_str());
  fprintf(f, "/>\n");

}