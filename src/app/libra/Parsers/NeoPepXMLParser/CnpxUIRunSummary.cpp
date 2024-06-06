#include "CnpxUIRunSummary.h"

using namespace std;


CnpxUIRunSummary::CnpxUIRunSummary(){
  pipelineIndex=0;
  runs=NULL;
}

CnpxUIRunSummary::~CnpxUIRunSummary(){
  runs=NULL;
}

CnpxMSMSRunSummary& CnpxUIRunSummary::operator[](const size_t& index){
  if (runs == NULL){
    cerr << "ERROR: CnpxUIRunSummary object is pointing to NULL." << endl;
    exit(-91);
  } else if (index >= runs->size()){
    cerr << "ERROR: Requested msms_run_summary beyond CnpxUIRunSummary boundary." << endl;
    exit(-92);
  }
  return runs->at(index);
}

size_t CnpxUIRunSummary::getPipelineIndex(){
  return pipelineIndex;
}

void CnpxUIRunSummary::set(std::vector<CnpxMSMSRunSummary>* p, size_t index){
  pipelineIndex=index;
  runs=p;
}

size_t CnpxUIRunSummary::size(){
  if (runs == NULL) return 0;
  else return runs->size();
}
