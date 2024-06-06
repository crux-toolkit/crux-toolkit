#ifndef _CNPXSEARCHSCORE_H
#define _CNPXSEARCHSCORE_H

#include "NeoPepXMLStructs.h"
#include <string>

class CnpxSearchScore {
public:

  void write(FILE* f, int tabs = -1);

  std::string name;
  std::string value;
 

private:

};

#endif
