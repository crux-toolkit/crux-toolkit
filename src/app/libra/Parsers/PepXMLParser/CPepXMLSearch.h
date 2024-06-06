#ifndef _CPEPXMLSEARCH_H
#define _CPEPXMLSEARCH_H

#include "PepXMLStructs.h"
#include <string>
#include <vector>

using namespace std;

class CPepXMLSearch {
public:

  CPepXMLSearch();
  CPepXMLSearch(const CPepXMLSearch& c);
  ~CPepXMLSearch();

  CPepXMLSearch& operator=(const CPepXMLSearch& c);

  void clear();
  size_t sizeOf();

  char DBindex;   //database searched
  
  string alg;     //algorithm used
  string version; //algorithm version
  
  vector<PepXMLSearchMod>*  mods;
  vector<PepXMLParam>*      params;

};

#endif
