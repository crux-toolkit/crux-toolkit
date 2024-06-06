#ifndef _CNPRINDISTINGUISHABLEPROTEIN_H
#define _CNPRINDISTINGUISHABLEPROTEIN_H

#include "NeoProtXMLStructs.h"
#include "CnprAnnotation.h"
#include "CnprParameter.h"
#include <string>
#include <vector>

class CnprIndistinguishableProtein {
public:

  void write(FILE* f, int tabs = -1);

  std::string protein_name;

  std::vector<CnprAnnotation> annotation;
  std::vector<CnprParameter> parameter;


private:

};

#endif