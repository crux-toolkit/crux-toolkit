#ifndef _CNPXISOTOPICCONTRIBUTIONS_H
#define _CNPXISOTOPICCONTRIBUTIONS_H

#include "NeoPepXMLStructs.h"
#include "CnpxContributingChannel.h"
#include <string>
#include <vector>

class CnpxIsotopicContributions {
public:

  void write(FILE* f, int tabs = -1);

  std::vector<CnpxContributingChannel> contributing_channel;

private:

};

#endif