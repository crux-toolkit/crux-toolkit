#ifndef PROT_ENZ_H
#define PROT_ENZ_H

#include <iostream>
#include <stdio.h>

#include "Common/constants.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Parsers/Algorithm2XML/SearchParams/SearchParams.h"
#include "Common/Array.h"
#include "Parsers/Parser/Tag.h"

class ProteolyticEnzyme {

 public:
  ProteolyticEnzyme(const char* name, const char* fidelity, Boolean indep, const char* description);
  ProteolyticEnzyme(Tag* tag);
  ~ProteolyticEnzyme();
  Boolean changeDefaultMinSpacing(int default_min_spacing);
  void enterSpecificity(const char* cut, const char* no_cut, Boolean n_sens, int min_spacing=1);
  void enterSpecificity(Tag* tag);
  void fixSpecificity();
  Array<Tag*>* getPepXMLTags();
  void write(ostream& os);
  void writePepXMLTags(ostream& os);
  void writeTraditionalPepXMLTags(FILE * fp);
  void writeTraditionalPepXMLTags(ogzstream * fp);
  int getNumTolTerm(char prev, const char* pep, char next);
  int getNumMissedCleavages(const char* pep, ModificationInfo* modinfo = NULL);
  Boolean hasSpecificity(char aa);
  Boolean hasModAASpecificity(SearchParams* params);
  char* strip(char* pep, Boolean remove_mods);
  const char* getName() const;
  const char* getFidelity() const;
  Boolean isSemiSpecific() const;
  int getMinNumberTermini() const;
  static char* getStandardEnzymeName(const char* input_name);
  static void parseEnzymeName(const char* input_name, char **name, Boolean &isSemiSpecific);


 protected:
  char* name_;
  char* fidelity_;
  Boolean independent_;
  char* description_;
  Array<char*>* cuts_;
  Array<char*>* no_cuts_;
  Array<Boolean>* n_sens_;
  Array<int>* min_spacing_;
  Boolean specified_;
  static char *default_fidelity_;
  static char default_description_[];
  static int default_min_spacing_;
  static Boolean default_independent_;
  static char default_nocut_[];

 public:
  static char CONSTANT_NAME_NONSPECIFIC_[];
  static char CONSTANT_FIDELITY_SEMISPECIFIC_[];
  static char CONSTANT_FIDELITY_SPECIFIC_[];
};

#endif
