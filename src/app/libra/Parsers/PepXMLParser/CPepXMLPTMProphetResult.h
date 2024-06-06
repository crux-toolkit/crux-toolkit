#ifndef _CPEPXMLPTMPROPHETRESULT_H
#define _CPEPXMLPTMPROPHETRESULT_H

#include "PepXMLStructs.h"
#include <string>
#include <vector>

class CPepXMLPTMProphetResult {
public:

  CPepXMLPTMProphetResult();
  CPepXMLPTMProphetResult(const CPepXMLPTMProphetResult& c);
  ~CPepXMLPTMProphetResult();

  CPepXMLPTMProphetResult& operator=(const CPepXMLPTMProphetResult& c);

  void clear();

  double prior;
  std::string ptm;
  std::string ptm_peptide;
  std::vector<PepXMLPepScore>* parameters;
  std::vector<PepXMLPTMMod>* mods;

};

#endif