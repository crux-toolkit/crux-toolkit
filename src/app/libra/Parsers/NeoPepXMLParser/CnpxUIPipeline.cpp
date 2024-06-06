#include "CnpxUIPipeline.h"

using namespace std;

CnpxUIPipeline::CnpxUIPipeline(){
  pipeline=NULL;
}

CnpxUIPipeline::~CnpxUIPipeline(){
  pipeline=NULL;
}

CnpxMSMSPipelineAnalysis& CnpxUIPipeline::operator[](const size_t& index){
  if (pipeline == NULL){
    cerr << "ERROR: CnpxUIPipeline object is pointing to NULL." << endl;
    exit(-81);
  } else if (index >= pipeline->size()){
    cerr << "ERROR: Requested msms_pipeline_analysis beyond CnpxUIPipeline boundary." << endl;
    exit(-82);
  }
  return pipeline->at(index);
}

void CnpxUIPipeline::set(std::vector<CnpxMSMSPipelineAnalysis>* p){
  pipeline=p;
}

size_t CnpxUIPipeline::size(){
  if (pipeline == NULL) return 0;
  else return pipeline->size();
}