#ifndef _CNPXMSMSPIPELINEANALYSIS_H
#define _CNPXMSMSPIPELINEANALYSIS_H

#include "CnpxAnalysisSummary.h"
#include "CnpxMSMSRunSummary.h"
#include "NeoPepXMLStructs.h"

#include <string>
#include <vector>

class CnpxMSMSPipelineAnalysis {
public:

  CnpxMSMSRunSummary* addMSMSRunSummary(std::string baseName, std::string rawDataType, std::string rawData);
  void write(FILE* f, int tabs=-1);

  std::vector<CnpxAnalysisSummary> analysis_summary;
  std::vector<CnpxMSMSRunSummary>  msms_run_summary;
  npxDateTime date;
  //CnpxDatasetDerivation dataset_derivation;
  std::string name;
  std::string summary_xml;
  

private:
  
};

#endif