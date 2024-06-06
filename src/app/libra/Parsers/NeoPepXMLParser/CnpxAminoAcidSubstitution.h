#ifndef _CNPXAMINOACIDSUBSTITUTION_H
#define _CNPXAMINOACIDSUBSTITUTION_H

#include "NeoPepXMLStructs.h"
#include <string>
#include <vector>

class CnpxAminoAcidSubstitution {
public:
  CnpxAminoAcidSubstitution();

  void write(FILE* f, int tabs = -1);

  int position;
  std::string orig_aa;
  int num_tol_term;
  std::string peptide_prev_aa;
  std::string peptide_next_aa;

private:

};

#endif
