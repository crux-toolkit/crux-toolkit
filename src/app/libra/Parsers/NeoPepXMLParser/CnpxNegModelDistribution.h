#ifndef _CNPXNEGMODELDISTRIBUTION_H
#define _CNPXNEGMODELDISTRIBUTION_H

#include "CnpxParameter.h"
#include <iostream>
#include <string>
#include <vector>

class CnpxNegModelDistribution {
public:
  CnpxNegModelDistribution();

  void write(FILE* f);

  std::string type;

  std::vector<CnpxParameter> parameter;

private:

};

#endif
