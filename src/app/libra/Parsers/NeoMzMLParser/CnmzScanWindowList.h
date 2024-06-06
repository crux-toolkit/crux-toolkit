#ifndef _CNMZSCANWINDOWLIST_H
#define _CNMZSCANWINDOWLIST_H

#include "NeoMzMLStructs.h"
#include "CnmzScanWindow.h"

#include <string>
#include <vector>

class CnmzScanWindowList {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::vector<CnmzScanWindow> scanWindow;

  int count;


private:

};

#endif