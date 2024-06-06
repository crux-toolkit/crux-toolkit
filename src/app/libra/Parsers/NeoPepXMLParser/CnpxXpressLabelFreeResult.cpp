#include "CnpxXpressLabelFreeResult.h"

using namespace std;

CnpxXpressLabelFreeResult::CnpxXpressLabelFreeResult() {
  charge=0;
  first_scan=-1;
  last_scan=-1;
  peak_intensity_scan=-1;
  active = false;
}

CnpxXpressLabelFreeResult::CnpxXpressLabelFreeResult(bool b) {
  charge = 0;
  first_scan = -1;
  last_scan = -1;
  peak_intensity_scan = -1;
  active = b;
}

bool CnpxXpressLabelFreeResult::present() {
  return active;
}

void CnpxXpressLabelFreeResult::write(FILE* f) {
  fprintf(f, "<xpresslabelfree_result charge=\"%d\"", charge);
  fprintf(f, " first_scan=\"%d\"",first_scan);
  fprintf(f, " last_scan=\"%d\"", last_scan);
  fprintf(f, " first_scan_RT_seconds=\"%s\"", first_scan_RT_seconds.c_str());
  fprintf(f, " last_scan_RT_seconds=\"%s\"", last_scan_RT_seconds.c_str());
  fprintf(f, " precursor_mz=\"%s\"", precursor_mz.c_str());
  fprintf(f, " peak_area=\"%s\"", peak_area.c_str());
  fprintf(f, " peak_intensity=\"%s\"", peak_intensity.c_str());
  fprintf(f, " peak_intensity_RT_seconds=\"%s\"", peak_intensity_RT_seconds.c_str());
  fprintf(f, " peak_intensity_scan=\"%d\"", peak_intensity_scan);
  fprintf(f,"/>\n");
}