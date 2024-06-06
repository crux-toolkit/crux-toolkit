#ifndef _CNPXSEARCHHIT_H
#define _CNPXSEARCHHIT_H

#include "CnpxAlternativeProtein.h"
#include "CnpxAnalysisResult.h"
#include "CnpxModificationInfo.h"
#include "CnpxSearchScore.h"
#include "CnpxXLink.h"
#include "NeoPepXMLStructs.h"
#include <string>
#include <vector>

class CnpxSearchHit {
public:
  CnpxSearchHit();

  CnpxAlternativeProtein* addAlternativeProtein(std::string protein);
  CnpxSearchScore* addSearchScore(std::string name, std::string value);
  std::string getModifiedPeptide();
  void write(FILE* f, int tabs=-1);

  int hit_rank;
  std::string peptide;
  std::string peptide_prev_aa;
  std::string peptide_next_aa;
  int peptide_start_pos;
  std::string protein;
  int protein_link_pos_a;
  int protein_link_pos_b;
  int num_tot_proteins;
  int num_matched_ions;
  int tot_num_ions;
  double calc_neutral_pep_mass;
  double massdiff;
  int num_tol_term;
  int num_missed_cleavages;
  int num_matched_peptides;
  std::string xlink_type;
  int is_rejected;
  std::string protein_descr;
  double calc_pI;
  double protein_mw;

  std::vector<CnpxModificationInfo> modification_info;  //0 or 1 object
  std::vector<CnpxXLink> xlink;                         //0 or 1 object

  std::vector<CnpxAlternativeProtein> alternative_protein;
  std::vector<CnpxAnalysisResult> analysis_result;
  std::vector<CnpxSearchScore> search_score;

private:

};

#endif
