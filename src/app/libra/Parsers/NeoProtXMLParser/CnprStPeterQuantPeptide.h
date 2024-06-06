#ifndef _CNPRSTPETERQUANTPEPTIDE_H
#define _CNPRSTPETERQUANTPEPTIDE_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprStPeterQuantPeptide {
public:

  CnprStPeterQuantPeptide();

  void write(FILE* f, int tabs = -1);

  std::string sequence;
  int charge;
  double dSI;
  double SI;
  double dSC;
  double SC;

private:

};

#endif