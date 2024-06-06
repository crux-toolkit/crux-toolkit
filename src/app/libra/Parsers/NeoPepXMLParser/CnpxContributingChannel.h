#ifndef _CNPXCONTRIBUTINGCHANNEL_H
#define _CNPXCONTRIBUTINGCHANNEL_H

#include "NeoPepXMLStructs.h"
#include "CnpxAffectedChannel.h"
#include <string>
#include <vector>

class CnpxContributingChannel {
public:
  CnpxContributingChannel();

  void write(FILE* f, int tabs = -1);

  int channel;

  std::vector<CnpxAffectedChannel> affected_channel;

private:

};

#endif