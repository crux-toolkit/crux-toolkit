#ifndef _CNPRCONTRIBUTINGCHANNEL_H
#define _CNPRCONTRIBUTINGCHANNEL_H

#include "NeoProtXMLStructs.h"
#include "CnprAffectedChannel.h"
#include <string>
#include <vector>

class CnprContributingChannel {
public:

  CnprContributingChannel();

  void write(FILE* f, int tabs = -1);

  int channel;
  std::vector<CnprAffectedChannel> affected_channel;

private:

};

#endif