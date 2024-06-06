#ifndef _CNPXXLINK_H
#define _CNPXXLINK_H

#include "CnpxLinkedPeptide.h"
#include "CnpxXLinkScore.h"
#include <string>

class CnpxXLink {
public:
  CnpxXLink();

  void write(FILE* f, int tabs = -1);

  std::string identifier;
  double mass;

  std::vector<CnpxLinkedPeptide> linked_peptide;
  std::vector<CnpxXLinkScore> xlink_score;

private:

};

#endif 
