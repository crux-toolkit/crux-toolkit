#ifndef _CPEPXMLANALYSISRESULT_H
#define _CPEPXMLANALYSISRESULT_H

#include "CPepXMLProphetResult.h"
#include "CPepXMLPTMProphetResult.h"
#include "PepXMLStructs.h"
#include <cstdio>
#include <string>
#include <vector>

enum eAnalysisResultType {
  empty = 0,
  PeptideProphet,
  InterProphet,
  PTMProphet
};

class CPepXMLAnalysisResult {
public:
  CPepXMLAnalysisResult();
  CPepXMLAnalysisResult(const CPepXMLAnalysisResult& c);
  ~CPepXMLAnalysisResult();

  CPepXMLAnalysisResult& operator=(const CPepXMLAnalysisResult& c);

  void del();
  void setAnalysisResult(void* ptr, eAnalysisResultType t);

  //for direct access to data
  eAnalysisResultType type;   //result structure type
  void* analysisResult;       //pointer to the analysis result structure/object
  std::vector<PepXMLPepScore>* parameters;  //the analysis result can have its own analysis parameters

private:

};

#endif