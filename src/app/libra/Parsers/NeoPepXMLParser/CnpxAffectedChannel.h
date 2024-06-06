#ifndef _CNPXAFFECTEDCHANNEL_H
#define _CNPXAFFECTEDCHANNEL_H

#include "NeoPepXMLStructs.h"
#include <string>
#include <vector>

class CnpxAffectedChannel {
public:
  CnpxAffectedChannel();

  void write(FILE* f, int tabs = -1);

  int channel;
  double correction;

private:

};

#endif