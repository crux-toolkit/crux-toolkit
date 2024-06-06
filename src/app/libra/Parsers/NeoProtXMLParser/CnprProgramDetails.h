#ifndef _CNPRPROGRAMDETAILS_H
#define _CNPRPROGRAMDETAILS_H

#include "NeoProtXMLStructs.h"
#include "CnprProteinProphetDetails.h"
#include <string>
#include <vector>

class CnprProgramDetails {
public:

  void write(FILE* f, int tabs = -1);

  std::string analysis;
  nprDateTime time;
  std::string version;

  std::vector<CnprProteinProphetDetails> proteinprophet_details;

private:

};

#endif