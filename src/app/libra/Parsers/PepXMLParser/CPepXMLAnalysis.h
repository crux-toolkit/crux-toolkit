#ifndef _CPEPXMLANALYSIS_H
#define _CPEPXMLANALYSIS_H

#include "PepXMLStructs.h"
#include <string>
#include <vector>

class CPepXMLAnalysis {
public:

  CPepXMLAnalysis();
  void clear();
  
  std::string analysis; //name of the analysis
  std::string version;

};

#endif
