#ifndef _CNPXQUANTICRESULT_H
#define _CNPXQUANTICRESULT_H

#include <string>
#include <vector>

class CnpxQuanticResult {
public:

  CnpxQuanticResult();
  CnpxQuanticResult(bool b);

  bool present();
  void write(FILE* f);

  double antic;

private:
  bool active;

};

#endif
