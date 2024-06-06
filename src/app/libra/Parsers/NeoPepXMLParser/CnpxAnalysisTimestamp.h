#ifndef _CNPXANALYSISTIMESTAMP_H
#define _CNPXANALYSISTIMESTAMP_H

#include "NeoPepXMLStructs.h"
#include "CnpxDatabaseRefreshTimestamp.h"
#include <iostream>
#include <string>

class CnpxAnalysisTimestamp {
public:
  CnpxAnalysisTimestamp();

  void write(FILE* f);

  std::string analysis;
  int id;
  npxDateTime time;

  CnpxDatabaseRefreshTimestamp database_refresh_timestamp;

private:

};

#endif 
