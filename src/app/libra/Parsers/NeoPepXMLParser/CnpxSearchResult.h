#ifndef _CNPXSEARCHRESULT_H
#define _CNPXSEARCHRESULT_H

#include "NeoPepXMLStructs.h"
#include "CnpxSearchHit.h"
#include <vector>

class CnpxSearchResult {
public:
  CnpxSearchResult();

  CnpxSearchHit* addSearchHit(int hitRank, std::string peptide, std::string protein, int numTotProteins, double calcNeutPepMass, double massdiff);
  void write(FILE* f, int tabs=-1);

  int search_id;
  
  std::vector<CnpxSearchHit> search_hit;

private:

};

#endif
