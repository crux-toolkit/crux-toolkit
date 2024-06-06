#ifndef _CNPXALTERNATIVEPROTEIN_H
#define _CNPXALTERNATIVEPROTEIN_H

#include "NeoPepXMLStructs.h"
#include <string>

class CnpxAlternativeProtein {
public:
  CnpxAlternativeProtein();

  void write(FILE* f, int tabs=-1);

  std::string protein;
  std::string protein_descr;
  int num_tol_term;
  double protein_mw;
  std::string peptide_prev_aa;
  std::string peptide_next_aa;
  int peptide_start_pos;
  int protein_link_pos_a;
  int protein_link_pos_b;

private:

};

#endif 
