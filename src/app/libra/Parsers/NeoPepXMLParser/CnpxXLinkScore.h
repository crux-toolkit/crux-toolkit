#ifndef _CNPXXLINKSCORE_H
#define _CNPXXLINKSCORE_H

#include "NeoPepXMLStructs.h"
#include <string>
#include <stdio.h>

class CnpxXLinkScore {
public:

  void write(FILE* f, int tabs=-1);

  std::string name;
  std::string type;
  std::string value;

private:

};

#endif 
