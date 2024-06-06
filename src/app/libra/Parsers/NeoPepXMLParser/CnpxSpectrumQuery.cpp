#include "CnpxSpectrumQuery.h"

using namespace std;

CnpxSpectrumQuery::CnpxSpectrumQuery(){
  spectrum.clear();
  spectrumNativeID.clear();
  start_scan=0;
  end_scan=0;
  retention_time_sec=0;
  collision_energy=0;
  compensation_voltage=0;
  precursor_intensity=0;
  activation_method.clear();
  precursor_neutral_mass=0;
  assumed_charge=0;
  search_specification.clear();
  index=0;
}

CnpxSearchResult* CnpxSpectrumQuery::addSearchResult(){
  CnpxSearchResult sr;
  search_result.push_back(sr);
  return &search_result.back();
}

void CnpxSpectrumQuery::write(FILE* f, int tabs){
  size_t i;

  string el = "spectrum_query";
  if (spectrum.empty()) NPXerrMsg(el, "spectrum");
  if (start_scan == 0) NPXerrMsg(el, "start_scan");
  if (end_scan == 0) NPXerrMsg(el, "end_scan");
  if (precursor_neutral_mass == 0) NPXerrMsg(el, "precursor_neutral_mass");
  if (assumed_charge == 0) NPXerrMsg(el, "assumed_charge");
  if (index == 0) NPXerrMsg(el, "index");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<spectrum_query spectrum=\"%s\"", spectrum.c_str());
  if(spectrumNativeID.size()>0) fprintf(f, " spectrumNativeID=\"%s\"", spectrumNativeID.c_str());
  fprintf(f, " start_scan=\"%d\"", start_scan);
  fprintf(f, " end_scan=\"%d\"", end_scan);
  fprintf(f, " precursor_neutral_mass=\"%.6lf\"", precursor_neutral_mass);
  fprintf(f, " assumed_charge=\"%d\"", assumed_charge);
  fprintf(f, " index=\"%d\"", index);
  if(retention_time_sec>0) fprintf(f, " retention_time_sec=\"%.1lf\"", retention_time_sec);
  if (collision_energy>0) fprintf(f, " collison_energy=\"%.1lf\"", collision_energy);
  if (compensation_voltage>0) fprintf(f, " compensation_voltage=\"%.2lf\"", compensation_voltage);
  if (precursor_intensity>0) fprintf(f, " precursor_intensity=\"%.2lf\"", precursor_intensity);
  if (!activation_method.empty()) fprintf(f, " activation_method=\"%s\"", activation_method.c_str());
  if (!search_specification.empty()) fprintf(f, " search_specification=\"%s\"", search_specification.c_str());
  fprintf(f,">\n");

  for(i=0;i<search_result.size();i++) search_result[i].write(f,t);
  
  NPXprintTabs(f, tabs);
  fprintf(f, "</spectrum_query>\n");

}
