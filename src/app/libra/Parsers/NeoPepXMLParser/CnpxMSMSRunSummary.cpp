#include "CnpxMSMSRunSummary.h"

using namespace std;


CnpxSearchSummary* CnpxMSMSRunSummary::addSearchSummary(std::string baseName, std::string searchEngine, std::string precursorMassType, std::string fragmentMassType, int searchID){
  CnpxSearchSummary s;
  s.base_name=baseName;
  s.search_engine=searchEngine;
  s.precursor_mass_type=precursorMassType;
  s.fragment_mass_type=fragmentMassType;
  s.search_id=searchID;
  search_summary.push_back(s);
  return &search_summary.back();
}

CnpxSpectrumQuery* CnpxMSMSRunSummary::addSpectrumQuery(std::string spec, int startScan, int endScan, double precursorNeutMass, int assumedCharge, int index){
  CnpxSpectrumQuery sq;
  sq.spectrum=spec;
  sq.start_scan=startScan;
  sq.end_scan=endScan;
  sq.precursor_neutral_mass=precursorNeutMass;
  sq.assumed_charge=assumedCharge;
  sq.index=index;
  spectrum_query.push_back(sq);
  return &spectrum_query.back();
}

void CnpxMSMSRunSummary::clear(){
  base_name.clear();
  msDetector.clear();
  msIonization.clear();
  msManufacturer.clear();
  msMassAnalyzer.clear();
  msModel.clear();
  raw_data.clear();
  raw_data_type.clear();
  sample_enzyme.clear();
  search_summary.clear();
  analysis_timestamp.clear();
  spectrum_query.clear();
}

void CnpxMSMSRunSummary::write(FILE* f, int tabs){
  size_t i;

  string el = "msms_run_summary";
  if (base_name.empty()) NPXerrMsg(el, "base_name");
  if (raw_data_type.empty()) NPXerrMsg(el, "raw_data_type");
  if (raw_data.empty()) NPXerrMsg(el, "raw_data");

  int t = tabs;
  if (t>-1) t++;
  NPXprintTabs(f, tabs);

  fprintf(f, "<msms_run_summary base_name=\"%s\"", base_name.c_str());
  fprintf(f, " raw_data_type=\"%s\"", raw_data_type.c_str());
  fprintf(f, " raw_data=\"%s\"", raw_data.c_str());
  if (!msManufacturer.empty()) fprintf(f, " msManufacturer=\"%s\"", msManufacturer.c_str());
  if (!msModel.empty()) fprintf(f, " msModel=\"%s\"", msModel.c_str());
  if (!msIonization.empty()) fprintf(f, " msIonization=\"%s\"", msIonization.c_str());
  if (!msMassAnalyzer.empty()) fprintf(f, " msMassAnalyzer=\"%s\"", msMassAnalyzer.c_str());
  if (!msDetector.empty()) fprintf(f, " msDetector=\"%s\"", msDetector.c_str());
  fprintf(f, ">\n");

  for (i = 0; i<sample_enzyme.size(); i++) sample_enzyme[i].write(f,t);
  for (i = 0; i<search_summary.size(); i++) search_summary[i].write(f,t);
  for(i=0;i<analysis_timestamp.size(); i++) analysis_timestamp[i].write(f);
  for(i=0;i<spectrum_query.size();i++) spectrum_query[i].write(f,t);

  NPXprintTabs(f, tabs);
  fprintf(f, "</msms_run_summary>\n");
}