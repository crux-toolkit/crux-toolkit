#include "CHardklorVariant.h"

using namespace std;

CHardklorVariant::CHardklorVariant(){
  ID=0;
  atoms = new vector<sInt>;
  enrich = new vector<sEnrichMercury>;
};

CHardklorVariant::CHardklorVariant(const CHardklorVariant& c){
  size_t i;
   
  atoms = new vector<sInt>;
  enrich = new vector<sEnrichMercury>;
  
  for(i=0;i<c.atoms->size();i++){
    atoms->push_back(c.atoms->at(i));
  };
  
  for(i=0;i<c.enrich->size();i++){
    enrich->push_back(c.enrich->at(i));
  };

  ID = c.ID;

};

CHardklorVariant::~CHardklorVariant(){
  delete atoms;
  delete enrich;
};

CHardklorVariant& CHardklorVariant::operator=(const CHardklorVariant& c){
  size_t i;
  if(this!=&c){
    delete atoms;
    delete enrich;
    atoms = new vector<sInt>;
    enrich = new vector<sEnrichMercury>;

    for(i=0;i<c.atoms->size();i++){
      atoms->push_back(c.atoms->at(i));
    };

    for(i=0;i<c.enrich->size();i++){
      enrich->push_back(c.enrich->at(i));
    };

    ID = c.ID;

  };
  return *this;
};

void CHardklorVariant::addAtom(const sInt& a){
  atoms->push_back(a);
};

void CHardklorVariant::addAtom(const int& a, const int& b){
  sInt s;
  s.iLower = a;
  s.iUpper = b;
  atoms->push_back(s);
};

void CHardklorVariant::addEnrich(const sEnrichMercury& a){
  enrich->push_back(a);
};

void CHardklorVariant::addEnrich(const int& a, const int& c, const double& b){
  sEnrichMercury s;
  s.atomNum = a;
  s.isotope = c;
  s.ape = b;
  enrich->push_back(s);
};

sInt& CHardklorVariant::atAtom(int i){
  return atoms->at(i);
};

sEnrichMercury& CHardklorVariant::atEnrich(int i){
  return enrich->at(i);
};

void CHardklorVariant::clear(){
  delete atoms;
  delete enrich;
  atoms = new vector<sInt>;
  enrich = new vector<sEnrichMercury>;
};

int CHardklorVariant::sizeAtom(){
  return atoms->size();
};

int CHardklorVariant::sizeEnrich(){
  return enrich->size();
};
