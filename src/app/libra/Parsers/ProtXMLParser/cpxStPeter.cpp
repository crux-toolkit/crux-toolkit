#include "cpxStPeter.h"

using namespace std;

//Constructors & Destructors
cpxStPeter::cpxStPeter(){
  counts=0;
  SI = 0;
  SIn = 0;
  NSAF = 0;
  ng = 0;
  ngC = 0;
  stpPeptide = new vector<spxSTPQPeptide>;
}

cpxStPeter::cpxStPeter(const cpxStPeter& c){
  counts = c.counts;
  SI = c.SI;
  SIn = c.SIn;
  NSAF = c.NSAF;
  ng = c.ng;
  ngC = c.ngC;
  stpPeptide = new vector<spxSTPQPeptide>;
  for(size_t i=0;i<c.stpPeptide->size();i++) stpPeptide->push_back(c.stpPeptide->at(i));
}

cpxStPeter::~cpxStPeter(){
  delete stpPeptide;
}

//Operators
cpxStPeter& cpxStPeter::operator=(const cpxStPeter& c){
  if(this!=&c){
    counts = c.counts;
    SI = c.SI;
    SIn = c.SIn;
    NSAF = c.NSAF;
    ng = c.ng;
    ngC = c.ngC;
    delete stpPeptide;
    stpPeptide = new vector<spxSTPQPeptide>;
    for (size_t i = 0; i<c.stpPeptide->size(); i++) stpPeptide->push_back(c.stpPeptide->at(i));
  }
  return *this;
}

spxSTPQPeptide& cpxStPeter::operator[](const size_t& index){
  return stpPeptide->at(index);
}
