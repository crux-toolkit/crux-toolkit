#include "CPepXMLSearch.h"

CPepXMLSearch::CPepXMLSearch(){
  DBindex=0;
  alg.clear();
  version.clear();
  mods=new vector<PepXMLSearchMod>;
  params=new vector<PepXMLParam>;
}

CPepXMLSearch::CPepXMLSearch(const CPepXMLSearch& c){
  //size_t i;
  DBindex = c.DBindex;
  alg = c.alg;
  version = c.version;
  mods = new vector<PepXMLSearchMod>(*c.mods);
  //for (i = 0; i<c.mods->size(); i++) mods->push_back(c.mods->at(i));
  params = new vector<PepXMLParam>(*c.params);
  //for (i = 0; i<c.params->size(); i++) params->push_back(c.params->at(i));
}

CPepXMLSearch::~CPepXMLSearch(){
  delete mods;
  delete params;
}

CPepXMLSearch& CPepXMLSearch::operator=(const CPepXMLSearch& c){
  if (this != &c){
    //size_t i;
    DBindex = c.DBindex;
    alg = c.alg;
    version = c.version;
    delete mods;
    delete params;
    mods = new vector<PepXMLSearchMod>(*c.mods);
    //for (i = 0; i<c.mods->size(); i++) mods->push_back(c.mods->at(i));
    params = new vector<PepXMLParam>(*c.params);
    //for (i = 0; i<c.params->size(); i++) params->push_back(c.params->at(i));
  }
  return *this;
}

void CPepXMLSearch::clear(){
  DBindex = 0;
  alg.clear();
  version.clear();
  delete mods;
  mods = new vector<PepXMLSearchMod>;
  delete params;
  params = new vector<PepXMLParam>;
}

size_t CPepXMLSearch::sizeOf(){
  size_t i;
  i=sizeof(*this);
  i+=alg.capacity();
  i+=version.capacity();
  i+=mods->size()*sizeof(PepXMLSearchMod);
  i += params->size()*sizeof(PepXMLParam);
  return i;
}