#ifndef _CNMZCV_H
#define _CNMZCV_H

#include "NeoMzMLStructs.h"
#include <string>


class CnmzCv {
public:

  void write(FILE* f, int tabs = -1, bool iterative = false);

  std::string URI;
  std::string fullName;
  std::string id;
  std::string version;

private:

};

#endif