#include "CPepXMLAnalysisResult.h"

using namespace std;

CPepXMLAnalysisResult::CPepXMLAnalysisResult(){
  type=empty;
  analysisResult=NULL;
  parameters=new vector<PepXMLPepScore>;
}

CPepXMLAnalysisResult::CPepXMLAnalysisResult(const CPepXMLAnalysisResult& c){
  CPepXMLProphetResult* pr;
  CPepXMLPTMProphetResult* ppr;

  type=c.type;
  switch(type){
  case PeptideProphet:
  case InterProphet:
    pr = new CPepXMLProphetResult(*(CPepXMLProphetResult*)c.analysisResult);
    analysisResult=pr;
    break;
  case PTMProphet:
    ppr = new CPepXMLPTMProphetResult(*(CPepXMLPTMProphetResult*)c.analysisResult);
    analysisResult = ppr;
    break;
  case empty:
  default:
    analysisResult=NULL;
    break;
  }

  parameters = new vector<PepXMLPepScore>(*c.parameters);
}

CPepXMLAnalysisResult::~CPepXMLAnalysisResult(){
  delete parameters;
  del();
}

CPepXMLAnalysisResult& CPepXMLAnalysisResult::operator=(const CPepXMLAnalysisResult& c){
  if(this!=&c){
    del();
    delete parameters;
    type = c.type;
    parameters = new vector<PepXMLPepScore>(*c.parameters);

    CPepXMLProphetResult* pr;
    CPepXMLPTMProphetResult* ppr;

    type = c.type;
    switch (type){
    case PeptideProphet:
    case InterProphet:
      pr = new CPepXMLProphetResult(*(CPepXMLProphetResult*)c.analysisResult);
      analysisResult = pr;
      break;
    case PTMProphet:
      ppr = new CPepXMLPTMProphetResult(*(CPepXMLPTMProphetResult*)c.analysisResult);
      analysisResult = ppr;
      break;
    case empty:
    default:
      analysisResult = NULL;
      break;
    }
  }
  return *this;
}

void CPepXMLAnalysisResult::setAnalysisResult(void* ptr, eAnalysisResultType t){
  //if current result exists, free it now; note that parameters remain the same.
  del();
  analysisResult=ptr;
  type=t;
}

void CPepXMLAnalysisResult::del(){
  if (analysisResult != NULL) {
    switch(type){
    case PeptideProphet:
    case InterProphet:
      delete (CPepXMLProphetResult*)analysisResult;
      break;
    case PTMProphet:
      delete (CPepXMLPTMProphetResult*)analysisResult;
      break;
    case empty:
      printf("Error: non-null analysis result\n");
      break;
    default:
      printf("Error: unknown analysis result\n");
      break;
    }
  }
  analysisResult=NULL;
  type=empty;
}
