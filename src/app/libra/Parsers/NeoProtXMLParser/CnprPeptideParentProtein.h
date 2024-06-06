#ifndef _CNPRPEPTIDEPARENTPROTEIN_H
#define _CNPRPEPTIDEPARENTPROTEIN_H

#include "NeoProtXMLStructs.h"
#include <string>
#include <vector>

class CnprPeptideParentProtein {
public:

  void write(FILE* f, int tabs = -1);

  std::string protein_name;

private:

};

#endif