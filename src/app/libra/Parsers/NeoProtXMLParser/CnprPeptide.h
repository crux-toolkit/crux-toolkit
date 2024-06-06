#ifndef _CNPRPEPTIDE_H
#define _CNPRPEPTIDE_H

#include "NeoProtXMLStructs.h"
#include "CnprIndistinguishablePeptide.h"
#include "CnprModificationInfo.h"
#include "CnprParameter.h"
#include "CnprPeptideParentProtein.h"
#include <string>
#include <vector>

class CnprPeptide {
public:

  CnprPeptide();

  void write(FILE* f, int tabs = -1);

  std::string peptide_sequence;
  int charge;
  double initial_probability;
  double nsp_adjusted_probability;
  double fpkm_adjusted_probability;
  double ni_adjusted_probability;
  double exp_sibling_ion_instances;
  double exp_sibling_ion_bin;
  double exp_tot_instances;
  std::string peptide_group_designator;
  double weight;
  std::string is_nondegenerate_evidence;
  int n_enzymatic_termini;
  double n_sibling_peptides;
  int n_sibling_peptides_bin;
  double max_fpkm;
  int fpkm_bin;
  int n_instances;
  double calc_neutral_pep_mass;
  std::string is_contributing_evidence;

  std::vector<CnprParameter> parameter;
  std::vector<CnprModificationInfo> modification_info;
  std::vector<CnprPeptideParentProtein> peptide_parent_protein;
  std::vector<CnprIndistinguishablePeptide> indistinguishable_peptide;


private:

};

#endif