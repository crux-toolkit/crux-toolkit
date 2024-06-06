#ifndef _CNPRPANALYSISSUMMARY_H
#define _CNPRPANALYSISSUMMARY_H

#include "NeoProtXMLStructs.h"
#include "CnprDecoyAnalysis.h"
#include "CnprDecoyAnalysisSummary.h"
#include "CnprLibraSummary.h"
#include "CnprStPeterAnalysisSummary.h"
#include <string>
#include <vector>

class CnprAnalysisSummary {
public:

  CnprAnalysisSummary();

  void write(FILE* f, int tabs = -1);

  std::string analysis;
  nprDateTime time;
  int id;

  std::vector<CnprDecoyAnalysisSummary> decoy_analysis_summary;
  std::vector<CnprDecoyAnalysis> decoy_analysis;
  std::vector<CnprLibraSummary> libra_summary;
  std::vector<CnprStPeterAnalysisSummary> StPeter_analysis_summary;

private:

};

#endif