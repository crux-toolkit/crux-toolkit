#ifndef _CNPRANNOTATION_H
#define _CNPRANNOTATION_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprAnnotation {
public:

  void write(FILE* f, int tabs = -1);

  std::string protein_description;
  std::string ipi_name;
  std::string refseq_name;
  std::string swissprot_name;
  std::string ensembl_name;
  std::string trembl_name;
  std::string locus_link_name;
  std::string flybase;


private:

};

#endif