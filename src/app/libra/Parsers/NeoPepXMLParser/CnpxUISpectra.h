#ifndef _CNPXUISPECTRA_H
#define _CNPXUISPECTRA_H

#include "CnpxSpectrumQuery.h"
#include <vector>

class CnpxUISpectra {
public:
  CnpxUISpectra();
  ~CnpxUISpectra();
  
  CnpxSpectrumQuery& operator[](const size_t& index);

  CnpxSearchHit& getHit(const size_t& queryIndex, const size_t& rank);
  size_t getPipelineIndex();
  size_t getRunSummaryIndex();
  void set(std::vector<CnpxSpectrumQuery>* p, const size_t pipeIndex, const size_t runIndex);
  size_t size();

private:
  size_t pipelineIndex;
  size_t runSummaryIndex;
  std::vector<CnpxSpectrumQuery>* spectra;

};

#endif 
