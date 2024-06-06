#ifndef _CPEPXMLPROPHETRESULT_H
#define _CPEPXMLPROPHETRESULT_H

#include "PepXMLStructs.h"
#include <string>
#include <vector>

class CPepXMLProphetResult {
public:

  CPepXMLProphetResult();
  CPepXMLProphetResult(const CPepXMLProphetResult& c);
  ~CPepXMLProphetResult();

  CPepXMLProphetResult& operator=(const CPepXMLProphetResult& c);

  void clear();

  double probability;
  double nttProbability[3];
  std::vector<PepXMLPepScore>* parameters;

};

#endif