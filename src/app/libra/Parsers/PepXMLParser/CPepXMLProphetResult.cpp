#include "CPepXMLProphetResult.h"

using namespace std;

CPepXMLProphetResult::CPepXMLProphetResult(){
  probability=0;
  nttProbability[0]=0;
  nttProbability[1] = 0;
  nttProbability[2] = 0;
  parameters = new vector<PepXMLPepScore>;
}

CPepXMLProphetResult::CPepXMLProphetResult(const CPepXMLProphetResult& c){
  probability = c.probability;
  nttProbability[0] = c.nttProbability[0];
  nttProbability[1] = c.nttProbability[1];
  nttProbability[2] = c.nttProbability[2];
  parameters = new vector<PepXMLPepScore>(*c.parameters);
}

CPepXMLProphetResult::~CPepXMLProphetResult(){
  delete parameters;
}

CPepXMLProphetResult& CPepXMLProphetResult::operator=(const CPepXMLProphetResult& c){
  if(this!=&c){
    probability = c.probability;
    nttProbability[0] = c.nttProbability[0];
    nttProbability[1] = c.nttProbability[1];
    nttProbability[2] = c.nttProbability[2];
    delete parameters;
    parameters = new vector<PepXMLPepScore>(*c.parameters);
  }
  return *this;
}

void CPepXMLProphetResult::clear(){
  probability = 0;
  nttProbability[0] = 0;
  nttProbability[1] = 0;
  nttProbability[2] = 0;
  delete parameters;
  parameters = new vector<PepXMLPepScore>;
}