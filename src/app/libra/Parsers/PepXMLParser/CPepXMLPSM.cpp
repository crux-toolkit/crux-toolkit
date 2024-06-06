#include "CPepXMLPSM.h"  

using namespace std;

CPepXMLPSM::CPepXMLPSM(){
  xlIndex=0;
  xlType=0;
  rank=0;
  calcPSMNeutMass=0;
  probability=0;
  iProphetProbability=0;
  peptide=NULL;
  xlPeptide=NULL;
  psmScores = new vector<PepXMLPepScore>;
  analysisResults = new vector<CPepXMLAnalysisResult>;
}

CPepXMLPSM::CPepXMLPSM(const CPepXMLPSM& p){
  xlIndex=p.xlIndex;
  xlType=p.xlType;
  rank=p.rank;
  calcPSMNeutMass=p.calcPSMNeutMass;
  probability=p.probability;
  iProphetProbability=p.iProphetProbability;
  peptide=NULL;
  if(p.peptide!=NULL) peptide = new CPepXMLPeptide(*p.peptide);
  xlPeptide=NULL;
  if(p.xlPeptide!=NULL) xlPeptide= new CPepXMLPeptide(*p.xlPeptide);
  psmScores = new vector<PepXMLPepScore>(*p.psmScores);
  analysisResults = new vector<CPepXMLAnalysisResult>(*p.analysisResults);
}

CPepXMLPSM::~CPepXMLPSM(){
  delete psmScores;
  delete analysisResults;
  if(peptide!=NULL) delete peptide;
  if(xlPeptide!=NULL) delete xlPeptide;
}

CPepXMLPSM& CPepXMLPSM::operator=(const CPepXMLPSM& p){
  if(this!=&p){
    if(peptide!=NULL) delete peptide;
    if(xlPeptide!=NULL) delete xlPeptide;
    delete psmScores;
    delete analysisResults;

    xlIndex=p.xlIndex;
    xlType=p.xlType;
    rank=p.rank;
    calcPSMNeutMass=p.calcPSMNeutMass;
    probability=p.probability;
    iProphetProbability=p.iProphetProbability;
    peptide=NULL;
    xlPeptide=NULL;
    if(p.peptide!=NULL) peptide = new CPepXMLPeptide(*p.peptide);
    if(p.xlPeptide!=NULL) xlPeptide= new CPepXMLPeptide(*p.xlPeptide);

    psmScores = new vector<PepXMLPepScore>(*p.psmScores);
    analysisResults = new vector<CPepXMLAnalysisResult>(*p.analysisResults);
  }
  return *this;
}

void CPepXMLPSM::clear(){
  xlIndex=0;
  xlType=0;
  rank=0;
  calcPSMNeutMass=0;
  probability=0;
  iProphetProbability=0;
  if(peptide!=NULL) delete peptide;
  if(xlPeptide!=NULL) delete xlPeptide;
  peptide=NULL;
  xlPeptide=NULL;
  delete psmScores;
  delete analysisResults;
  psmScores = new vector<PepXMLPepScore>;
  analysisResults = new vector<CPepXMLAnalysisResult>;
}

size_t CPepXMLPSM::sizeOf(){
  size_t i=sizeof(*this)+psmScores->capacity()*sizeof(PepXMLPepScore);
  i+=analysisResults->capacity()*sizeof(CPepXMLAnalysisResult);
  if(peptide!=NULL) i+=peptide->sizeOf();
  if(xlPeptide!=NULL) i+=xlPeptide->sizeOf();
  return i;
}