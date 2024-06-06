#ifndef _CNPXMSMSRUNSUMMARY_H
#define _CNPXMSMSRUNSUMMARY_H

#include "CnpxAnalysisTimestamp.h"
#include "CnpxSampleEnzyme.h"
#include "CnpxSearchSummary.h"
#include "CnpxSpectrumQuery.h"
#include <string>
#include <vector>

class CnpxMSMSRunSummary {
public:

  CnpxSearchSummary* addSearchSummary(std::string baseName, std::string searchEngine, std::string precursorMassType, std::string fragmentMassType, int searchID);
  CnpxSpectrumQuery* addSpectrumQuery(std::string spec, int startScan, int endScan, double precursorNeutMass, int assumedCharge, int index);
  void clear();
  void write(FILE* f, int tabs=-1);

  std::string base_name;
  std::string msDetector;
  std::string msIonization;
  std::string msManufacturer;
  std::string msMassAnalyzer;
  std::string msModel;
  std::string raw_data;
  std::string raw_data_type;

  std::vector<CnpxSampleEnzyme> sample_enzyme;
  std::vector<CnpxSearchSummary> search_summary;
  std::vector<CnpxAnalysisTimestamp> analysis_timestamp;
  std::vector<CnpxSpectrumQuery> spectrum_query;

private:

};

#endif
