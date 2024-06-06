#ifndef _CNPXSEARCHSUMMARY_H
#define _CNPXSEARCHSUMMARY_H

#include "CnpxAminoAcidModification.h"
#include "CnpxEnzymaticSearchConstraint.h"
#include "CnpxParameter.h"
#include "CnpxSearchDatabase.h"
#include "CnpxSequenceSearchConstraint.h"
#include "CnpxTerminalModification.h"
#include <string>
#include <vector>


class CnpxSearchSummary {
public:

  CnpxSearchSummary();

  CnpxAminoAcidModification* addAminoAcidModification(std::string aminoAcid, double massDiff, double mass, std::string variable);
  CnpxEnzymaticSearchConstraint* addEnzymaticSearchConstraint(std::string enzyme, int maxInternalCleavages, int minTermini);
  CnpxParameter* addParameter(std::string name, std::string value);
  CnpxSearchDatabase* addSearchDatabase(std::string localPath, std::string type);
  void write(FILE* f, int tabs=-1);

  std::string base_name;
  std::string fragment_mass_type;
  std::string out_data;
  std::string out_data_type;
  std::string precursor_mass_type;
  std::string search_engine;
  std::string search_engine_version;
  int search_id;

  std::vector<CnpxAminoAcidModification> aminoacid_modification;
  std::vector<CnpxEnzymaticSearchConstraint> enzymatic_search_constraint;
  std::vector<CnpxParameter> parameter;
  std::vector<CnpxSearchDatabase> search_database;
  std::vector<CnpxSequenceSearchConstraint> sequence_search_constraint;
  std::vector<CnpxTerminalModification> terminal_modification;

private:

};

#endif