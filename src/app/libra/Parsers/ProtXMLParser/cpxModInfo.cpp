#include "cpxModInfo.h"

using namespace std;

cpxModInfo::cpxModInfo(){
  modifiedPeptide.clear();
  modAminoAcidMass = new vector<spxMod>;
  nTermMod=0;
  cTermMod=0;
}

cpxModInfo::cpxModInfo(const cpxModInfo& c){
  modifiedPeptide=c.modifiedPeptide;
  nTermMod=c.nTermMod;
  cTermMod=c.cTermMod;
  modAminoAcidMass = new vector<spxMod>;
  for (size_t i = 0; i<c.modAminoAcidMass->size(); i++) modAminoAcidMass->push_back(c.modAminoAcidMass->at(i));
}

cpxModInfo::~cpxModInfo(){
  delete modAminoAcidMass;
}

cpxModInfo& cpxModInfo::operator=(const cpxModInfo& c){
  if (this != &c){
    modifiedPeptide = c.modifiedPeptide;
    nTermMod = c.nTermMod;
    cTermMod = c.cTermMod;
    delete modAminoAcidMass;
    modAminoAcidMass = new vector<spxMod>;
    for (size_t i = 0; i<c.modAminoAcidMass->size(); i++) modAminoAcidMass->push_back(c.modAminoAcidMass->at(i));
  }
  return *this;
}
