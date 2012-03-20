#include "CHardklorProtein.h"
#include <string.h>

using namespace std;

CHardklorProtein::CHardklorProtein(){
  enrich = new vector<sEnrichMercury>;
}

CHardklorProtein::CHardklorProtein(const CHardklorProtein& c){
  size_t i;

  enrich = new vector<sEnrichMercury>;
  for(i=0;i<c.enrich->size();i++) enrich->push_back(c.enrich->at(i));

  ID = c.ID;
  mz = c.mz;
  monoMass = c.monoMass;
  shft = c.shft;
  abun = c.abun;
  charge = c.charge;
  C = c.C;
  rTime = c.rTime;
  strcpy(seq,c.seq);

}

CHardklorProtein::~CHardklorProtein(){
  delete enrich;
}

CHardklorProtein& CHardklorProtein::operator=(const CHardklorProtein& c){
  size_t i;

  if(this != &c){
    delete enrich;
    enrich = new vector<sEnrichMercury>;
    for(i=0;i<c.enrich->size();i++) enrich->push_back(c.enrich->at(i));
    
    ID = c.ID;
    mz = c.mz;
    monoMass = c.monoMass;
    shft = c.shft;
    abun = c.abun;
    charge = c.charge;
    C = c.C;
    rTime = c.rTime;
    strcpy(seq,c.seq);
  }
  return *this;

}

void CHardklorProtein::add(int a, int c, double b){
  sEnrichMercury s;
  s.atomNum = a;
  s.isotope = c;
  s.ape = b;
  enrich->push_back(s);
}

sEnrichMercury& CHardklorProtein::at(int i){
  return enrich->at(i);
}

void CHardklorProtein::clear(){
  delete enrich;
  enrich = new vector<sEnrichMercury>;
}

int CHardklorProtein::size(){
  return enrich->size();
}
