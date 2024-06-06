#ifndef _CNPXENZYMATICSEARCHCONSTRAINT_H
#define _CNPXENZYMATICSEARCHCONSTRAINT_H

#include <iostream>
#include <string>

class CnpxEnzymaticSearchConstraint {
public:
  CnpxEnzymaticSearchConstraint();

  void write(FILE* f);

  std::string enzyme;
  int max_num_internal_cleavages;
  int min_number_termini;

private:

};

#endif