#ifndef _CNPXLINKEDPEPTIDE_H
#define _CNPXLINKEDPEPTIDE_H

#include "CnpxAlternativeProtein.h"
#include "CnpxModificationInfo.h"
#include "CnpxXLinkScore.h"
#include <string>

class CnpxLinkedPeptide {
public:
  CnpxLinkedPeptide();

  void write(FILE* f, int tabs=-1);

  double calc_neutral_pep_mass;
  double complement_mass;
  std::string designation;
  int num_tot_proteins;
  std::string peptide;
  std::string peptide_prev_aa;
  std::string peptide_next_aa;
  int peptide_start_pos;
  std::string protein;
  int protein_link_pos_a;

  std::vector<CnpxAlternativeProtein> alternative_protein;
  std::vector<CnpxModificationInfo> modification_info;
  std::vector<CnpxXLinkScore> xlink_score;

private:

};

#endif 
