#ifndef _CNPXUIPIPELINE_H
#define _CNPXUIPIPELINE_H

#include "CnpxMSMSPipelineAnalysis.h"
#include <vector>

class CnpxUIPipeline {
public:
  CnpxUIPipeline();
  ~CnpxUIPipeline();

  CnpxMSMSPipelineAnalysis& operator[](const size_t& index);

  void set(std::vector<CnpxMSMSPipelineAnalysis>* p);
  size_t size();

private:
  std::vector<CnpxMSMSPipelineAnalysis>* pipeline;

};

#endif 
