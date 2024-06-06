#ifndef _CNPRSTPETERQUANT_H
#define _CNPRSTPETERQUANT_H

#include "NeoProtXMLStructs.h"
#include "CnprStPeterQuantPeptide.h"
#include <string>
#include <vector>

class CnprStPeterQuant {
public:

  CnprStPeterQuant();

  void write(FILE* f, int tabs = -1);

  double dSI;
  double dSIn;
  double SI;
  double SIn;
  double dCounts;
  double counts;
  double dNSAF;
  double NSAF;
  double ng;
  double ngC;

  std::vector<CnprStPeterQuantPeptide> StPeterQuant_peptide;

private:

};

#endif