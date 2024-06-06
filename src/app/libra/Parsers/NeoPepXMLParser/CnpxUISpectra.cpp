#include "CnpxUISpectra.h"

using namespace std;

CnpxUISpectra::CnpxUISpectra(){
  pipelineIndex=0;
  runSummaryIndex=0;
  spectra = NULL;
}

CnpxUISpectra::~CnpxUISpectra(){
  spectra=NULL;
}

CnpxSpectrumQuery& CnpxUISpectra::operator[](const size_t& index){
  if(spectra==NULL){
    cerr << "ERROR, CnpxUISpectra::operator[]: CnpxUISpectra object is pointing to NULL." << endl;
    exit(-101);
  } else if(index>=spectra->size()){
    cerr << "ERROR, CnpxUISpectra::operator[]: Requested spectrum beyond CnpxUISpectra boundary." << endl;
    exit(-101);
  }
  return spectra->at(index);
}

CnpxSearchHit& CnpxUISpectra::getHit(const size_t& queryIndex, const size_t& rank){
  if(rank<1){
    cerr << "ERROR, CnpxUISpectra::getHit(): user requested rank less than 1." << endl;
    exit(-101);
  }
  if(queryIndex>=spectra->size()){
    cerr << "ERROR, CnpxUISpectra::getHit(): spectrum_query index out of bounds." << endl;
    exit(-101);
  } else if(rank>spectra->at(queryIndex).search_result[0].search_hit.size()){
    cerr << "ERROR, CnpxUISpectra::getHit(): requested rank out of bounds." << endl;
    exit(-101);
  }
  return spectra->at(queryIndex).search_result[0].search_hit[rank-1];
}

size_t CnpxUISpectra::getPipelineIndex(){
  return pipelineIndex;
}

size_t CnpxUISpectra::getRunSummaryIndex(){
  return runSummaryIndex;
}

void CnpxUISpectra::set(std::vector<CnpxSpectrumQuery>* p, const size_t pipeIndex, const size_t runIndex){
  pipelineIndex=pipeIndex;
  runSummaryIndex=runIndex;
  spectra=p;
}

size_t CnpxUISpectra::size(){
  if(spectra==NULL) return 0;
  else return spectra->size();
}