#ifndef _CPEPXMLSPECTRUM_H
#define _CPEPXMLSPECTRUM_H

#include "CPepXMLPSM.h"

class CPepXMLSpectrum {
public:

  CPepXMLSpectrum();
  CPepXMLSpectrum(const CPepXMLSpectrum& p);
  ~CPepXMLSpectrum();

  CPepXMLSpectrum& operator=(const CPepXMLSpectrum& p);
  CPepXMLPSM& operator[ ](const size_t& index);

  void    clear();
  size_t  size();
  size_t  sizeOf();

  char searchID;

  int charge;
  int fileID;
  int scanEnd;
  int scanNumber;
  int scanStart;

  float rTimeSec;

  double precursorNeutMass;

  std::string ID;
  std::string nativeID;

  std::vector<CPepXMLPSM>* psms;
};

#endif