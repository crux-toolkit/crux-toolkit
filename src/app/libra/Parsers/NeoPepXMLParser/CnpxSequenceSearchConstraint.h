#ifndef _CNPXSEQUENCESEARCHCONSTRAINT_H
#define _CNPXSEQUENCESEARCHCONSTRAINT_H

#include <iostream>
#include <string>

class CnpxSequenceSearchConstraint {
public:
  CnpxSequenceSearchConstraint();

  void write(FILE* f);

  std::string sequence;

private:

};

#endif