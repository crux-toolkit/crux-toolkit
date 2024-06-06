#include "CPepXMLPTMProphetResult.h"

using namespace std;

CPepXMLPTMProphetResult::CPepXMLPTMProphetResult(){
  prior=0;
  ptm.clear();
  ptm_peptide.clear();
  parameters = new vector<PepXMLPepScore>;
  mods = new vector<PepXMLPTMMod>;
}

CPepXMLPTMProphetResult::CPepXMLPTMProphetResult(const CPepXMLPTMProphetResult& c){
  prior = c.prior;
  ptm = c.ptm;
  ptm_peptide = c.ptm_peptide;
  parameters = new vector<PepXMLPepScore>(*c.parameters);
  mods = new vector<PepXMLPTMMod>(*c.mods);
}

CPepXMLPTMProphetResult::~CPepXMLPTMProphetResult(){
  delete parameters;
  delete mods;
}

CPepXMLPTMProphetResult& CPepXMLPTMProphetResult::operator=(const CPepXMLPTMProphetResult& c){
  if(this!=&c){
    prior = c.prior;
    ptm = c.ptm;
    ptm_peptide = c.ptm_peptide;
    delete parameters;
    delete mods;
    parameters = new vector<PepXMLPepScore>(*c.parameters);
    mods = new vector<PepXMLPTMMod>(*c.mods);
  }
  return *this;
}

void CPepXMLPTMProphetResult::clear(){
  delete parameters;
  delete mods;

  prior = 0;
  ptm.clear();
  ptm_peptide.clear();
  parameters = new vector<PepXMLPepScore>;
  mods = new vector<PepXMLPTMMod>;
}

