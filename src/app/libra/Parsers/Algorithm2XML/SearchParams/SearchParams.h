#ifndef SEARCH_PARAMS_H
#define SEARCH_PARAMS_H

//#define WRITE_MOD_INFO // to write mod info

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "Common/constants.h"
#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "Parsers/Parser/Tag.h"
#include "Common/ModificationInfo/ModificationInfo.h"

class Modification
{
public:
   Modification() {
      memset(this,0,sizeof(*this));
   }
  char aa;
  double massdiff;
  Boolean variable;
  char symbol[4];
  double mass;
  Boolean terminal;
  Boolean protein_terminus;
};



class SearchParams {
  friend class CombineOut;
 public:

  SearchParams(const char* paramsfile);
  virtual ~SearchParams();
  
  void setMassType(Boolean monoisotopic);
  Tag* getModificationTag(Modification* mod);
  Array<Tag*>* getModificationTags();
  Tag* getEnzymeConstraintTag();
  Tag* getSequenceConstraintTag(const char* seq);
  Array<Tag*>* getSequenceConstraintTags();
  char getModifiedAA(int index);
  int getNumModifiedAAs();

  virtual Array<Tag*>* getParameterTags() = 0;
  Array<Tag*>* getSearchParamTagsPlus(const char* basename, const char* engine, Boolean parent_mono, Boolean fragment_mono, const char* database);
  Array<Tag*>* getSearchParamTags(const char* basename, const char* engine, const char* database);
  const char* getParamsEnzyme();
  double getAverageAAMass(char aa);
  double getMonoisotopicAAMass(char aa);

  void setModificationSymbol(Modification* mod, Boolean advance);
  virtual void modifyAminoacidModifications(char* peptide) = 0;
  virtual void modifyTerminalModifications(char* peptide) = 0;

  Array<Tag*>* getModificationInfoTags(const char* peptide);
  char* getStandardModifiedPeptide(const char* peptide);

  Boolean getPrecursorMonoisotopic() { return precursor_monoisotopic_; }

 protected:
  
  virtual void init() = 0;
  virtual const char* getEnzyme(int index) = 0;
  virtual void setModificationMasses() = 0;
  void replaceModificationSymbol(Modification* mod, Boolean advance);
  void updateModifications();

  void replace(char* peptide, const char* orig, const char* replacements, Boolean termini);
  void setModificationSymbols(int num_sequence_mods, const char** sequence_modification_symbols, int num_terminal_mods, const char** terminal_modification_symbols);

  Array<Modification*>* modifications_;
  char* paramsfile_;

  Array<char*>* sequence_constraints_;

  //char sequence_constraint_[25];
  int enzyme_index_;
  int max_num_internal_cls_;
  int min_num_tol_term_;
  Boolean precursor_monoisotopic_;
  Boolean fragment_monoisotopic_;

  char** sequence_modification_symbols_;
  int num_sequence_modifications_;
  int sequence_modification_index_;
  char** terminal_modification_symbols_;
  int num_terminal_modifications_;
  int terminal_modification_index_;

};



#endif
