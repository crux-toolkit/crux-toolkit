#include "cpxProtein.h"

using namespace std;

cpxProtein::cpxProtein(){
  confidence=0;
  groupSiblingID.clear();
  nIndistinguishableProteins=0;
  pctSpectrumIDs=0;
  percentCoverage=0;
  probability=0;
  proteinName.clear();
  totalNumberDistinctPeptides=0;
  totalNumberPeptides=0;
  peptide = new vector<cpxPeptide>;
  uniqueStrippedPeptides = new vector<string>;
}

cpxProtein::cpxProtein(const cpxProtein& c){
  confidence = c.confidence;
  groupSiblingID = c.groupSiblingID;
  nIndistinguishableProteins = c.nIndistinguishableProteins;
  pctSpectrumIDs = c.pctSpectrumIDs;
  percentCoverage = c.percentCoverage;
  probability = c.probability;
  proteinName = c.proteinName;
  stPeter = c.stPeter;
  totalNumberDistinctPeptides = c.totalNumberDistinctPeptides;
  totalNumberPeptides = c.totalNumberPeptides;
  peptide = new vector<cpxPeptide>;
  uniqueStrippedPeptides = new vector<string>;
  size_t i;
  for (i = 0; i<c.peptide->size(); i++) peptide->push_back(c.peptide->at(i));
  for (i = 0; i<c.uniqueStrippedPeptides->size(); i++) uniqueStrippedPeptides->push_back(c.uniqueStrippedPeptides->at(i));
}

cpxProtein::~cpxProtein(){
  delete peptide;
  delete uniqueStrippedPeptides;
}

cpxPeptide& cpxProtein::operator[](const size_t index){
  return peptide->at(index);
}

cpxProtein& cpxProtein::operator=(const cpxProtein& c){
  if (this != &c){
    confidence = c.confidence;
    groupSiblingID = c.groupSiblingID;
    nIndistinguishableProteins = c.nIndistinguishableProteins;
    pctSpectrumIDs = c.pctSpectrumIDs;
    percentCoverage = c.percentCoverage;
    probability = c.probability;
    proteinName = c.proteinName;
    stPeter = c.stPeter;
    totalNumberDistinctPeptides = c.totalNumberDistinctPeptides;
    totalNumberPeptides = c.totalNumberPeptides;
    delete peptide;
    delete uniqueStrippedPeptides;
    peptide = new vector<cpxPeptide>;
    uniqueStrippedPeptides = new vector<string>;
    size_t i;
    for (i = 0; i<c.peptide->size(); i++) peptide->push_back(c.peptide->at(i));
    for (i = 0; i<c.uniqueStrippedPeptides->size(); i++) uniqueStrippedPeptides->push_back(c.uniqueStrippedPeptides->at(i));
  }
  return *this;
}

size_t cpxProtein::size(){
  return peptide->size();
}

