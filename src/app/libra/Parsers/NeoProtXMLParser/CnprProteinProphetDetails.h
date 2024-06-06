#ifndef _CNPRPROTEINPROPHETDETAILS_H
#define _CNPRPROTEINPROPHETDETAILS_H

#include "NeoProtXMLStructs.h"
#include "CnprErrorPoint.h"
#include "CnprNSPInformation.h"
#include "CnprProteinSummaryDataFilter.h"
#include <string>
#include <vector>

class CnprProteinProphetDetails {
public:

  void write(FILE* f, int tabs = -1);

  std::string occam_flag;
  std::string groups_flag;
  std::string degen_flag;
  std::string nsp_flag;
  std::string fpkm_flag;
  std::string initial_peptide_wt_iters;
  std::string nsp_distribution_iters;
  std::string final_peptide_wt_iters;
  std::string run_options;

  CnprNSPInformation nsp_information;
  //std::vector<CnprFPKMInformation> fpkm_information;
  //std::vector<CnprNIInformation> ni_information;
  std::vector<CnprProteinSummaryDataFilter> protein_summary_data_filter;
  std::vector<CnprErrorPoint> error_point;

private:

};

#endif