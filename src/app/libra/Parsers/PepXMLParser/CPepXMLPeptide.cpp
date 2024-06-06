#include "CPepXMLPeptide.h"

using namespace std;

CPepXMLPeptide::CPepXMLPeptide(){
  label=0;
  nextAA='-';
  prevAA='-';
  calcPepNeutMass=0;
  complementMass=0;
  modifiedPeptide.clear();
  peptide.clear();
  modifiedPeptide.shrink_to_fit();
  peptide.shrink_to_fit();
  mods=NULL; //new vector<PepXMLPepMod>;
  proteins=new vector<size_t>;
  xlScores=new vector<PepXMLPepScore>;
}

CPepXMLPeptide::CPepXMLPeptide(const CPepXMLPeptide& p){
  label=p.label;
  nextAA=p.nextAA;
  prevAA=p.prevAA;
  calcPepNeutMass=p.calcPepNeutMass;
  complementMass=p.complementMass;
  modifiedPeptide=p.modifiedPeptide;
  peptide=p.peptide;
  modifiedPeptide.shrink_to_fit();
  peptide.shrink_to_fit();
  if(p.mods!=NULL) mods=new vector<PepXMLPepMod>(*p.mods);
  else mods=NULL;
  proteins=new vector<size_t>(*p.proteins);
  xlScores=new vector<PepXMLPepScore>(*p.xlScores);
}

CPepXMLPeptide::~CPepXMLPeptide(){
  if(mods!=NULL) delete mods;
  delete proteins;
  delete xlScores;
}

CPepXMLPeptide& CPepXMLPeptide::operator=(const CPepXMLPeptide& p){
  if(this!=&p){
    label=p.label;
    nextAA=p.nextAA;
    prevAA=p.prevAA;
    calcPepNeutMass=p.calcPepNeutMass;
    complementMass=p.complementMass;
    modifiedPeptide=p.modifiedPeptide;
    peptide=p.peptide;
    modifiedPeptide.shrink_to_fit();
    peptide.shrink_to_fit();
    if(mods!=NULL) delete mods;
    delete proteins;
    delete xlScores;
    if(p.mods!=NULL) mods=new vector<PepXMLPepMod>(*p.mods);
    else mods=NULL;
    proteins=new vector<size_t>(*p.proteins);
    xlScores=new vector<PepXMLPepScore>(*p.xlScores);
  }
  return *this;
}

void CPepXMLPeptide::addMod(PepXMLPepMod& m){
  if(mods==NULL) mods=new vector<PepXMLPepMod>;
  mods->push_back(m);
}

void CPepXMLPeptide::clear(){
  label=0;
  nextAA='-';
  prevAA='-';
  calcPepNeutMass=0;
  complementMass=0;
  modifiedPeptide.clear();
  peptide.clear();
  modifiedPeptide.shrink_to_fit();
  peptide.shrink_to_fit();
  if(mods!=NULL) delete mods;
  delete proteins;
  delete xlScores;
  mods=NULL;
  proteins=new vector<size_t>;
  xlScores=new vector<PepXMLPepScore>;
}

size_t CPepXMLPeptide::sizeMod(){
  if(mods==NULL) return 0;
  else return mods->size();
}

size_t CPepXMLPeptide::sizeOf(){
  size_t i;
  i=sizeof(*this);
  i+=modifiedPeptide.capacity();
  i+=peptide.capacity();
  if (mods != NULL) i += mods->capacity()*sizeof(PepXMLPepMod);
  i+=proteins->capacity()*sizeof(size_t);
  i += xlScores->capacity()*sizeof(PepXMLPepScore);
  return i;
}

