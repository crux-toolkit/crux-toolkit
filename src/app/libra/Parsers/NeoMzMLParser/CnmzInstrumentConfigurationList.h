#ifndef _CNMZINSTRUMENTCONFIGURATIONLIST_H
#define _CNMZINSTRUMENTCONFIGURATIONLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzInstrumentConfiguration.h"
#include <vector>


class CnmzInstrumentConfigurationList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzInstrumentConfiguration> instrumentConfiguration;

  int count;

private:

};

#endif