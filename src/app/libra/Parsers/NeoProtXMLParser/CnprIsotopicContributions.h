#ifndef _CNPRISOTOPICCONTRIBUTIONS_H
#define _CNPRISOTOPICCONTRIBUTIONS_H

#include "NeoProtXMLStructs.h"
#include "CnprContributingChannel.h"
#include <string>
#include <vector>

class CnprIsotopicContributions {
public:

  void write(FILE* f, int tabs = -1);

  std::vector<CnprContributingChannel> contributing_channel;

private:

};

#endif