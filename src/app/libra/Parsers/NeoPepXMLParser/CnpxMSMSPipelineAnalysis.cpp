#include "CnpxMSMSPipelineAnalysis.h"

using namespace std;


CnpxMSMSRunSummary* CnpxMSMSPipelineAnalysis::addMSMSRunSummary(std::string baseName, std::string rawDataType, std::string rawData){
  CnpxMSMSRunSummary rs;
  rs.base_name=baseName;
  rs.raw_data_type=rawDataType;
  rs.raw_data=rawData;
  msms_run_summary.push_back(rs);
  return &msms_run_summary.back();
}

void CnpxMSMSPipelineAnalysis::write(FILE* f, int tabs){
  string el = "msms_pipeline_analysis";
  if(date.date.year==0) NPXerrMsg(el, "date");
  if (summary_xml.empty()) NPXerrMsg(el, "summary_xml");

  size_t i;
  fprintf(f, "<msms_pipeline_analysis date=\"%4d-%02d-%02dT%02d:%02d:%02d\"", date.date.year, date.date.month, date.date.day, date.time.hour, date.time.minute, date.time.second);
  fprintf(f, " xmlns=\"%s\"",npx_xmlns.c_str());
  fprintf(f, " xmlns:xsi=\"%s\"", npx_xmlns_xsi.c_str());
  fprintf(f, " xsi:schemaLocation=\"%s\"", npx_xsi_schemaLocation.c_str());
  fprintf(f, " summary_xml=\"%s\">\n", summary_xml.c_str());

  int t = tabs;
  if(t>-1) t++;

  for(i=0;i<analysis_summary.size();i++) analysis_summary[i].write(f,t);
  for(i=0;i<msms_run_summary.size();i++) msms_run_summary[i].write(f,t);

  fprintf(f, "</msms_pipeline_analysis>\n");
}