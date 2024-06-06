#ifndef _CNPRPROTEIN_H
#define _CNPRPROTEIN_H

#include "NeoProtXMLStructs.h"
#include "CnprAnalysisResult.h"
#include "CnprAnnotation.h"
#include "CnprIndistinguishableProtein.h"
#include "CnprParameter.h"
#include "CnprPeptide.h"
#include <string>
#include <vector>

class CnprProtein {
public:

  CnprProtein();

  void write(FILE* f, int tabs = -1);

  double confidence;
  std::string group_sibling_id;
  int n_indistinguishable_proteins;
  std::string pct_spectrum_ids;
  double percent_coverage;
  double probability;
  std::string protein_name;
  std::string subsuming_protein_entry;
  int total_number_distinct_peptides;
  int total_number_peptides;
  std::string unique_stripped_peptides;
  

  std::vector<CnprAnalysisResult> analysis_result;
  std::vector<CnprAnnotation> annotation;
  std::vector<CnprIndistinguishableProtein> indistinguishable_protein;
  std::vector<CnprParameter> parameter;
  std::vector<CnprPeptide> peptide;


private:

};

#endif