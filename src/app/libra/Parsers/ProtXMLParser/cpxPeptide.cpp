#include "cpxPeptide.h"

using namespace std;

cpxPeptide::cpxPeptide(){
  charge = 0;
  expTotInstances=0;
  fpkmAdjustedProbability=0;
  initialProbability=0;
  isContributingEvidence=0;
  isNondegenerateEvidence=0;
  nEnzymaticTermini=0;
  nInstances=0;
  nspAdjustedProbability=0;
  nSiblingPeptides=0;
  nSiblingPeptidesBin=0;
  peptideSequence.clear();
  weight=0;
  indistinguishablePeptide = new vector<cpxIndPeptide>;
}

cpxPeptide::cpxPeptide(const cpxPeptide& c){
  charge = c.charge;
  expTotInstances = c.expTotInstances;
  fpkmAdjustedProbability = c.fpkmAdjustedProbability;
  initialProbability = c.initialProbability;
  isContributingEvidence = c.isContributingEvidence;
  isNondegenerateEvidence = c.isNondegenerateEvidence;
  modificationInfo=c.modificationInfo;
  nEnzymaticTermini = c.nEnzymaticTermini;
  nInstances = c.nInstances;
  nspAdjustedProbability = c.nspAdjustedProbability;
  nSiblingPeptides = c.nSiblingPeptides;
  nSiblingPeptidesBin = c.nSiblingPeptidesBin;
  peptideSequence = c.peptideSequence;
  weight = c.weight;
  indistinguishablePeptide = new vector<cpxIndPeptide>;
  for (size_t i = 0; i<c.indistinguishablePeptide->size(); i++) indistinguishablePeptide->push_back(c.indistinguishablePeptide->at(i));
}

cpxPeptide::~cpxPeptide(){
  delete indistinguishablePeptide;
}

cpxIndPeptide& cpxPeptide::operator[](const size_t index){
  return indistinguishablePeptide->at(index);
}

cpxPeptide& cpxPeptide::operator=(const cpxPeptide& c){
  if (this != &c){
    charge = c.charge;
    expTotInstances = c.expTotInstances;
    fpkmAdjustedProbability = c.fpkmAdjustedProbability;
    initialProbability = c.initialProbability;
    isContributingEvidence = c.isContributingEvidence;
    isNondegenerateEvidence = c.isNondegenerateEvidence;
    modificationInfo = c.modificationInfo;
    nEnzymaticTermini = c.nEnzymaticTermini;
    nInstances = c.nInstances;
    nspAdjustedProbability = c.nspAdjustedProbability;
    nSiblingPeptides = c.nSiblingPeptides;
    nSiblingPeptidesBin = c.nSiblingPeptidesBin;
    peptideSequence = c.peptideSequence;
    weight = c.weight;
    delete indistinguishablePeptide;
    indistinguishablePeptide = new vector<cpxIndPeptide>;
    for (size_t i = 0; i<c.indistinguishablePeptide->size(); i++) indistinguishablePeptide->push_back(c.indistinguishablePeptide->at(i));
  }
  return *this;
}

size_t cpxPeptide::size(){
  return indistinguishablePeptide->size();
}
