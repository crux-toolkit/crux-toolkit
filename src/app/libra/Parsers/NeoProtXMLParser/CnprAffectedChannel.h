#ifndef _CNPRAFFECTEDCHANNEL_H
#define _CNPRAFFECTEDCHANNEL_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprAffectedChannel {
public:

  CnprAffectedChannel();

  void write(FILE* f, int tabs = -1);

  int channel;
  double correction;

private:

};

#endif