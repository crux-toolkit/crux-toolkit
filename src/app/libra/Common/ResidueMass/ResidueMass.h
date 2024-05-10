#ifndef RESIDUE_H
#define RESIDUE_H

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include "Common/constants.h"

class ResidueMass {

 public:

  ResidueMass();

  static double getMass(char res, Boolean monoisotopic);
  static double getProteinMass(const char* prot, Boolean monoisotopic);
  static double getStdModMass(const char* mod, Boolean monoisotopic, char label);
  static const char* getStdModResidues(const char* mod);
  static const char* getStdModResidues(const char* mod, Boolean& n_term_aa_mod, Boolean& c_term_aa_mod);


 protected:

  static double monoisotopic_masses_[];
  static double average_masses_[];
  static double monoisotopic_nterm_mass_;
  static double average_nterm_mass_;
  static double monoisotopic_cterm_mass_;
  static double average_cterm_mass_;
  static double getAverageMass(char res);
  static double getMonoisotopicMass(char res);


};


#endif
