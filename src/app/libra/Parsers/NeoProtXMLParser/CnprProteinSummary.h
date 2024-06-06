#ifndef _CNPRPROTEINSUMMARY_H
#define _CNPRPROTEINSUMMARY_H

#include "NeoProtXMLStructs.h"
#include "CnprAnalysisSummary.h"
#include "CnprDatasetDerivation.h"
#include "CnprProteinGroup.h"
#include "CnprProteinSummaryHeader.h"
#include <string>
#include <vector>

class CnprProteinSummary {
public:

  void write(FILE* f, int tabs=-1);

  std::string xmlns;
  std::string xmlns_xsi;
  std::string xsi_schemaLocation;
  std::string summary_xml;

  CnprProteinSummaryHeader protein_summary_header;
  std::vector<CnprAnalysisSummary> analysis_summary;
  CnprDatasetDerivation dataset_derivation;
  std::vector<CnprProteinGroup> protein_group;

private:

};

#endif