#ifndef _CNMZPRECURSOR_H
#define _CNMZPRECURSOR_H

#include "NeoMzMLStructs.h"
#include "CnmzActivation.h"
#include "CnmzIsolationWindow.h"
#include "CnmzSelectedIonList.h"
#include <vector>
#include <string>


class CnmzPrecursor {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzIsolationWindow> isolationWindow;
  std::vector<CnmzSelectedIonList> selectedIonList;
  CnmzActivation activation;

  std::string externalSpectrumID;
  std::string sourceFileRef;
  std::string spectrumRef;

private:

};

#endif