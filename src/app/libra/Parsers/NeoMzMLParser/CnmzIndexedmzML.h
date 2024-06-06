#ifndef _CNMZINDEXEDMZML_H
#define _CNMZINDEXEDMZML_H

#include "NeoMzMLStructs.h"
#include "CnmzFileChecksum.h"
#include "CnmzIndexList.h"
#include "CnmzIndexListOffset.h"
#include "CnmzMzML.h"

#include <string>
#include <vector>

class CnmzIndexedmzML {
public:

  CnmzIndexedmzML(){mzML=NULL;}
  ~CnmzIndexedmzML(){mzML=NULL;}

  void write(FILE* f, int tabs = -1, bool iterative=false);

  CnmzMzML* mzML;
  CnmzIndexList indexList;
  CnmzIndexListOffset indexListOffset;
  CnmzFileChecksum fileChecksum;

  std::string xmlns;
  std::string xmlns_xsi;
  std::string xsi_schemaLocation;


private:

};

#endif