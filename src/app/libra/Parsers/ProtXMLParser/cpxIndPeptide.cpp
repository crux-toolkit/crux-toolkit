#include "cpxIndPeptide.h"

cpxIndPeptide::cpxIndPeptide(){
  charge=0;
  calcNeutralPepMass=0;
  peptideSequence.clear();
}

cpxIndPeptide::cpxIndPeptide(const cpxIndPeptide& c){
  charge = c.charge;
  calcNeutralPepMass = c.calcNeutralPepMass;
  peptideSequence = c.peptideSequence;
  modificationInfo = c.modificationInfo;
}

cpxIndPeptide::~cpxIndPeptide(){
}

cpxIndPeptide& cpxIndPeptide::operator=(const cpxIndPeptide& c){
  if (this != &c){
    charge = c.charge;
    calcNeutralPepMass = c.calcNeutralPepMass;
    peptideSequence = c.peptideSequence;
    modificationInfo = c.modificationInfo;
  }
  return *this;
}
